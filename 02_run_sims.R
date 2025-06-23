library(dplyr)
library(tidyr)
library(survival)
library(nimble)
library(posterior)
library(bayesplot)
library(cmdstanr)
library(mclust)

source('/proj/ibrahimlab/leap/sims_pwe/models.R')

grid = expand.grid(p = c(0.0, 0.5, 1.0),
                   n = c(25, 50, 100),
                   q = c(-0.4, -0.2, 0, 0.2, 0.4))

id=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

for (hist in 1:nrow(grid)) {
  if (!is.na(id) && id%%nrow(grid)!=hist%%nrow(grid)) next # each job works on 1 histdata
  for (cur in 1:10000) {
    if (!is.na(id) && !(ceiling(id/nrow(grid))*20-19<=cur & cur<=ceiling(id/nrow(grid))*20)) next
    # each job works on 20 curdata, so total MAKE SURE TOTAL #jobs is nrow(grid)*10000/20
    if (grid$p[hist] == 1) next # no need to run p=1 scenario
    print(paste0('working on histdata = ', hist, ', curdata = ', cur))
    
    histdata.all = readRDS(paste0('data/historical/', 'histdata_', hist, '.rds'))
    data.all = readRDS(paste0('data/current/', 'curdata_', grid$n[hist], '_', cur, '.rds'))
    
    S = max(data.all$strata)
    n <- n0 <- grid$n[hist] %>% rep(S)
    J = numeric(length(n)); J[n==25] = 4; J[n==50] = 7; J[n==100] = 10
    K = rep(2, S)
    
    setup = sleap_setup(
      formula = Surv(y, event) ~ a, data = data.all%>%select(-c(exch)), histdata = histdata.all%>%select(-c(exch,k)), stratavar = 'strata'
      , nclasses = rep(2, S)
      , nintervals = J, pooled_intervals = FALSE
      , bhm = T
    )
    
    setup2 = sleap_setup(
      formula = Surv(y, event) ~ a, data = data.all%>%select(-c(exch)), histdata = histdata.all%>%select(-c(exch,k)), stratavar = 'strata'
      , nclasses = rep(2, S)
      , nintervals = J, pooled_intervals = FALSE
      , bhm = F
    )
    
    #----------------------------------#
    # Data processing for other models
    #----------------------------------
    
    a = histdata.all %>% rename(s=strata, trt=a, death=event) %>% group_by(s) %>% mutate(s_i = row_number(), histdata=1)
    b = data.all %>% rename(s=strata, trt=a, death=event) %>% group_by(s) %>% mutate(s_i = row_number(), histdata=0)
    c = bind_rows(a, b)
    
    # breakpoints are stratum-specific and based on current data
    c.long = data.frame()
    for (s in 1:S) {
      nbreaks = J[s]
      if (nbreaks==1) {
        breaks = c(0, 10000)
      } else {
        breaks = c(0, quantile(b$y[b$death==1 & b$s==s], probs = (1:(nbreaks-1)) / nbreaks), 10000)
      }
      c.long_s = survSplit(Surv(y, death) ~ ., cut = breaks, data = c[c$s==s,], end = 'tstop') %>%
        mutate(risktime = tstop - tstart, intervalind = NA)
      for ( j in seq_along(breaks) ) {
        c.long_s$intervalind <- ifelse(c.long_s$tstart == breaks[j], j, c.long_s$intervalind)
      }
      c.long = rbind(c.long, c.long_s)
    }
    
    event.long0=array(rep(0, S*max(n0)*max(J)), dim=c(S,max(n0),max(J)))
    logrisktime.long0=array(rep(-Inf, S*max(n0)*max(J)), dim=c(S,max(n0),max(J)))
    event.long=array(rep(0, S*max(n)*max(J)), dim=c(S,max(n),max(J)))
    logrisktime.long=array(rep(-Inf, S*max(n)*max(J)), dim=c(S,max(n),max(J)))
    for (row in 1:nrow(c.long)) {
      if (c.long$histdata[row]==1) {
        event.long0[c.long$s[row], c.long$s_i[row], c.long$intervalind[row]] = c.long$death[row]
        logrisktime.long0[c.long$s[row], c.long$s_i[row], c.long$intervalind[row]] = log(c.long$risktime[row])
      } else {
        event.long[c.long$s[row], c.long$s_i[row], c.long$intervalind[row]] = c.long$death[row]
        logrisktime.long[c.long$s[row], c.long$s_i[row], c.long$intervalind[row]] = log(c.long$risktime[row])
      }
    }
    
    Ji0 = array(rep(0, S*max(n0)), dim=c(S,max(n0)))
    x0 = array(rep(0, S*max(n0)), dim=c(S,max(n0)))
    Ji = array(rep(0, S*max(n)), dim=c(S,max(n)))
    x = array(rep(0, S*max(n)), dim=c(S,max(n)))
    for (row in 1:nrow(c)) {
      if (c$histdata[row]==1) {
        Ji0[c$s[row], c$s_i[row]] = sum(c.long$s==c$s[row] & c.long$s_i==c$s_i[row] & c.long$histdata==c$histdata[row])
        x0[c$s[row], c$s_i[row]] = c$trt[row]
      } else {
        Ji[c$s[row], c$s_i[row]] = sum(c.long$s==c$s[row] & c.long$s_i==c$s_i[row] & c.long$histdata==c$histdata[row])
        x[c$s[row], c$s_i[row]] = c$trt[row]
      }
    }
    
    nimdata  <- list(
      event = event.long, logrisktime = logrisktime.long, x = x, Ji=Ji
      , event0 = event.long0, logrisktime0 = logrisktime.long0, x0 = x0, Ji0=Ji0
    )
    
    #--------------------------------#
    # Run models
    #--------------------------------
    
    for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
      fname = paste0('results/raw/', method, '/res_', hist, '_', cur, '.rds')
      if (file.exists(fname) & file.size(fname)!=0) next
      
      if (method=='sleap') {
        #------------------#
        # SLEAP
        #------------------
        model <- nimbleModel(
          sleap, constants = setup$const, data = setup$data, inits = setup$inits
        )
        cmodel_sleap <- compileNimble(model)
        mcmcConf <- configureMCMC(model, print = FALSE, monitors = c('beta', 'beta_mean', 'beta_sd', 'blhaz', 'gamma'))
        mcmcConf$removeSampler(c('beta', 'beta_sd'))
        for ( s in 1:setup$const$S ) mcmcConf$addSampler( paste0('beta[', s, ']'), 'slice')
        mcmcConf$addSampler('beta_sd', 'slice')
        mcmcConf
        mcmc_sleap <- buildMCMC(mcmcConf)
        cmcmc    <- compileNimble(mcmc_sleap, project = cmodel_sleap, resetFunctions = TRUE)
        smpl     <- runMCMC(cmcmc, 2000 + 1 * 10000, 2000, thin = 1)
        rm(model, cmodel_sleap, mcmcConf, mcmc_sleap, cmcmc)
      } else if (method=='pp') {
        #------------------#
        # POWER PRIOR
        #------------------
        niminit  <- list(
          beta = setup$inits$beta
          , blhaz = setup$inits$blhaz[,1,]
          , beta_mean = setup$inits$beta_mean
          , beta_sd   = setup$inits$beta_sd
          , a0 = rep(1,S)
        )
        nimconst <- list(
          S=S, J=J, n=n, n0=n0
          , beta_mean_mean=setup$data$beta_mean_mean, beta_mean_sd=setup$data$beta_mean_sd
          , beta_sd_mu = setup$data$beta_sd_mu, beta_sd_sigma = setup$data$beta_sd_sigma
          , blhaz_shape=setup$data$blhaz_shape, blhaz_rate=setup$data$blhaz_rate
        )
        model    <- nimbleModel(pp, nimconst, nimdata, niminit)
        cmodel_pp   <- compileNimble(model)
        mcmcConf <- configureMCMC(model)
        mcmcConf$removeSampler('beta')
        for ( s in 1:setup$const$S ) mcmcConf$addSampler( paste0('beta[', s, ']'), 'slice')
        mcmcConf$removeSampler('beta_sd')
        mcmcConf$addSampler('beta_sd', 'slice')
        mcmcConf$addMonitors('beta', 'blhaz')
        mcmc_pp     <- buildMCMC(mcmcConf)
        cmcmc    <- compileNimble(mcmc_pp, project = cmodel_pp, resetFunctions = TRUE)
        smpl     <- runMCMC(cmcmc, 2000 + 1 * 10000, 2000, thin = 1)
        rm(model, cmodel_pp, mcmcConf, mcmc_pp, cmcmc)
      } else if (method=='bhm') {
        #------------------#
        # BHM
        #------------------
        niminit  <- list(
          beta = setup$inits$beta
          , blhaz = setup$inits$blhaz[,1,]
          , beta_mean = setup$inits$beta_mean
          , beta_sd   = setup$inits$beta_sd
        )
        nimconst <- list(
          S=S, J=J, n=n, n0=n0
          , beta_mean_mean=setup$data$beta_mean_mean, beta_mean_sd=setup$data$beta_mean_sd
          , beta_sd_mu = setup$data$beta_sd_mu, beta_sd_sigma = setup$data$beta_sd_sigma
          , blhaz_shape=setup$data$blhaz_shape, blhaz_rate=setup$data$blhaz_rate
        )
        model    <- nimbleModel(bhm, nimconst, nimdata, niminit)
        cmodel_bhm   <- compileNimble(model)
        mcmcConf <- configureMCMC(model)
        mcmcConf$removeSampler('beta')
        for ( s in 1:setup$const$S ) mcmcConf$addSampler( paste0('beta[', s, ']'), 'slice')
        mcmcConf$removeSampler('beta_sd')
        mcmcConf$addSampler('beta_sd', 'slice')
        mcmcConf$addMonitors('beta')
        mcmc_bhm     <- buildMCMC(mcmcConf)
        cmcmc    <- compileNimble(mcmc_bhm, project = cmodel_bhm, resetFunctions = TRUE)
        smpl     <- runMCMC(cmcmc, 2000 + 1 * 10000, 2000, thin = 1)
        rm(model, cmodel_bhm, mcmcConf, mcmc_bhm, cmcmc)
      } else if (method=='noninf') {
        #------------------#
        # NON-INFORMATIVE
        #------------------
        niminit  <- list(
          beta = setup$inits$beta
          , blhaz = setup$inits$blhaz[,1,]
        )
        nimconst <- list(
          S=S, K=K, J=J, n=n
          , beta_mean = 0, beta_sd   = 5
          , blhaz_shape=setup$data$blhaz_shape, blhaz_rate=setup$data$blhaz_rate
        )
        model    <- nimbleModel(noninf, nimconst, nimdata, niminit)
        cmodel_noninf   <- compileNimble(model)
        mcmcConf <- configureMCMC(model)
        mcmcConf$removeSampler('beta')
        for ( s in 1:setup$const$S ) mcmcConf$addSampler( paste0('beta[', s, ']'), 'slice')
        mcmcConf$addMonitors('beta')
        mcmc_noninf     <- buildMCMC(mcmcConf)
        cmcmc    <- compileNimble(mcmc_noninf, project = cmodel_noninf, resetFunctions = TRUE)
        smpl     <- runMCMC(cmcmc, 2000 + 1 * 10000, 2000, thin = 1)
        rm(model, cmodel_noninf, mcmcConf, mcmc_noninf, cmcmc)
      } else if (method=='sleap2') {
        #------------------#
        # SLEAP2 (no BHM)
        #------------------
        model <- nimbleModel(
          sleap2, constants = setup2$const, data = setup2$data, inits = setup2$inits
        )
        cmodel_sleap2 <- compileNimble(model)
        mcmcConf <- configureMCMC(model, print = FALSE, monitors = c('beta', 'blhaz', 'gamma'))
        mcmcConf$removeSampler(c('beta'))
        for ( s in 1:setup$const$S ) mcmcConf$addSampler( paste0('beta[', s, ']'), 'slice')
        mcmc_sleap2 <- buildMCMC(mcmcConf)
        cmcmc    <- compileNimble(mcmc_sleap2, project = cmodel_sleap2, resetFunctions = TRUE)
        smpl     <- runMCMC(cmcmc, 2000 + 1 * 10000, 2000, thin = 1)
        rm(model, cmodel_sleap2, mcmcConf, mcmc_sleap2, cmcmc)
      }
      
      
      if (method=='sleap2') {
        beta_mean = smpl[,sapply(1:S, function(s) paste0('beta[', s, ']'))] %>% rowMeans
        beta_sd = smpl[,sapply(1:S, function(s) paste0('beta[', s, ']'))] %>% apply(1, sd)
        smpl = cbind(beta_mean, beta_sd, smpl)
      }
      
      smpl.mcmc <- as_draws_matrix(smpl)
      # smpl.mcmc %>% mcmc_trace('beta[1]')
      # smpl.mcmc %>% mcmc_hist('beta[1]')
      res = 
        smpl.mcmc %>% summarize_draws(
          'mean', 'median', 'sd', ~quantile(.x, probs = c(0.025, 0.975))
          , 'ess_bulk'
        ) %>%
        mutate(across(where(is.numeric), round, 3))
      
      # res %>% filter(variable %in% sapply(1:S, function(s) paste0('gamma[', s, ', 1]')))
  
      saveRDS(res, file=paste0('results/raw/', method, '/res_', hist, '_', cur, '.rds'))
      
    }
  }
}