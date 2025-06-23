library(dplyr)
library(tidyr)
library(survival)
library(nimble)
library(posterior)
library(bayesplot)
library(cmdstanr)
library(mclust)
library(ggplot2)
library(tables)
library(MatchIt)

id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
gene = ifelse(id%%2==0, 2, 1)
run.ps = id>2
nsample = 250000
runBHMonly = F
pooled = F

source('models.R')

a = read.csv('data/surv_sanofi.csv') %>% rename_all(tolower) %>% rename(subjid=rsubjid) %>%
  mutate(histdata=1, kras=NA, race=case_when(race=='CAUCASIAN/WHITE' ~ 'White or Caucasian',
                                             race=='BLACK' ~ 'Black',
                                             race=='OTHER' ~ 'Other'),
         tumor_colon = as.factor(ifelse(tumor=='COLON', 'Colon', 'Other')),
         race_white = as.factor(ifelse(race=='White or Caucasian', race, 'Other')))
b = read.csv('data/surv_amgen.csv') %>% rename_all(tolower) %>% mutate(histdata=0) %>%
  filter(kras==ifelse(gene==1, 'Wild-type', 'Mutant')) %>%
  # filter(pfs!=0) %>% # only use for PFS
  mutate(kras=ifelse(kras=='', 'NA', kras),
         tumor_colon = as.factor(ifelse(tumor=='Colon', 'Colon', 'Other')),
         race_white = as.factor(ifelse(race=='White or Caucasian', race, 'Other')))
c = bind_rows(a, b) %>%
  mutate(across(all_of(c('ecog', 'ecog2', 'bevacizumab', 'oxaliplatin', 'histdata', 'kras', 'race', 'sex')), as.factor),
         # trt = factor(trt, level=c('1', '0')), # comment out if not doing descriptive stats.
         s = as.numeric(interaction(ecog2, bevacizumab, oxaliplatin)),
         curdata = as.factor(1 - as.numeric(as.character(histdata))))
a = c %>% filter(histdata==1)
b = c %>% filter(histdata==0)
        ###NOTE STRATA NUMBER HERE WILL LATER BE SORTED BY HISTDATA SIZE
# interaction(c$ecog2, c$bevacizumab, c$oxaliplatin)['Levels'] # see strata numbers

# --------------------------------
# Trimming using Propensity Score
# --------------------------------

if (run.ps) {
  common_vars = c('age', 'sex', 'race_white', 'tumor_colon')
  
  c.complete = c %>% filter(if_all(all_of(common_vars), ~ !is.na(.)))
  
  ps.fit = glm(as.formula(paste0('curdata ~ ', paste(common_vars, collapse = ' + '))),
               data = c.complete,
               family = 'binomial')
  
  ps = predict(ps.fit, type = "response")
  
  # Trimming 
  ps.cur.min = min(ps[c.complete$curdata==1])
  ps.cur.max = max(ps[c.complete$curdata==1])
  
  c.trimmed = c.complete %>%
    cbind(ps) %>%
    filter(histdata==1 | (ps>ps.cur.min & ps<ps.cur.max))
  
  # Matching
  
  matching = 
    matchit(as.formula(paste0('curdata ~ ', paste(common_vars, collapse = ' + '))),
            data = c.trimmed,
            method = "nearest",
            distance = "logit")
  
  c.matched = match.data(matching)
  a = c.matched %>% filter(histdata==1)
  c = bind_rows(a%>%select(names(b)), b)
}

# -----------------------------

a = c %>% filter(histdata==1) %>% group_by(s) %>% mutate(s_i = row_number()) %>% ungroup
b = c %>% filter(histdata==0) %>% group_by(s) %>% mutate(s_i = row_number()) %>% ungroup
c = bind_rows(a, b)

# If needed, remove empty strata/collapse strata number to make sure theres no gap
if (gene==1) {
  a = a %>% mutate(s=ifelse(s>4, s-1, s))
  b = b %>% mutate(s=ifelse(s>4, s-1, s))
  c = c %>% mutate(s=ifelse(s>4, s-1, s))
}

if (pooled) b = rbind(b, a)

S=max(c$s)
K=rep(2, S)

j_dat = readRDS(paste0('data_analysis_findBestJ_', ifelse(gene==1, 'WT', 'MUTANT'), '.rds')) %>%
  as.data.frame() %>%
  mutate(V1 = ifelse(is.na(V1), 0, V1))
J = apply(j_dat, 1, which.min) %>% unlist
# if (gene==1) {
#   J = c(8,1,3,6,1,5,1)
# } else {
#   J = c(9,1,1,1,8,1,3,1)
# }

n0 = a %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)
n = b %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)
a %>% group_by(s) %>% summarize(events=sum(death), n0=n())
b %>% group_by(s) %>% summarize(events=sum(death), n=n())

setup = sleap_setup(
  formula = Surv(os, death) ~ trt, data = b%>%rename(strata='s'), histdata = a%>%rename(strata='s'), stratavar = 'strata'
  , nclasses = rep(2, S)
  , nintervals = J, pooled_intervals = FALSE
  , bhm = T
)


#----------------------------------#
# Data processing for other models
#----------------------------------

# breakpoints are stratum-specific and based on current data
c.long = data.frame()
for (s in 1:S) {
  nbreaks = J[s]
  if (nbreaks==1) {
    breaks = c(0, 10000)
  } else {
    breaks = c(0, quantile(b$os[b$death==1 & b$s==s], probs = (1:(nbreaks-1)) / nbreaks), 10000)
  }
  c.long_s = survSplit(Surv(os, death) ~ ., cut = breaks, data = c[c$s==s,], end = 'tstop') %>%
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

res = list()
smpls = list()

for (method in c('sleap', 'pp', 'bhm', 'noninf', 'cwpp')) {
  if (runBHMonly & method!='rmapp.50') next
  if (pooled & method!='bhm') next
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
    smpl     <- runMCMC(cmcmc, 2000 + 1 * nsample, 2000, thin = 1)
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
    smpl     <- runMCMC(cmcmc, 2000 + 1 * nsample, 2000, thin = 1)
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
    smpl     <- runMCMC(cmcmc, 2000 + 1 * nsample, 2000, thin = 1)
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
    smpl     <- runMCMC(cmcmc, 2000 + 1 * nsample, 2000, thin = 1)
    rm(model, cmodel_noninf, mcmcConf, mcmc_noninf, cmcmc)
  } else if (grepl('^rmapp', method)) {
    #------------------#
    # ROBUST MAP PRIOR
    #------------------
    
    if (method=='rmapp.25') w=0.25
    if (method=='rmapp.50') w=0.50
    if (method=='rmapp.75') w=0.75
    
    if (runBHMonly==T) {
      # 1. BHM portion
      model_bhm = cmdstan_model('bhm_pwe.stan')
      smpl_bhm <- model_bhm$sample(
        data = list(
          'S' = S, 'J' = J, 'n' = n0, 'Ji' = Ji0,
          'event' = event.long0, 'logrisktime' = logrisktime.long0,
          'logblhaz_mean_mean' = 0, 'logblhaz_mean_sd' = 10,
          'logblhaz_sd_mean' = 0, 'logblhaz_sd_sd' = 1
        )
        , iter_warmup = 2000, iter_sampling = 12000, chains = 1, parallel_chains = 1
      )
      smpl.df <- posterior::as_draws_df(smpl_bhm$draws()) %>% as.data.frame()
      # colMeans(smpl.df)%>% head(30) %>% exp()
      # hist(exp(smpl.df$`logblhaz_mean[1]`))
      
      # 2. Fitting normal mixture using `Mclust`
      logblhaz = array(rep(0, nrow(smpl.df)*S*max(J)), dim=c(nrow(smpl.df),S,max(J)))
      G_max = 5
      G = rep(0, S)
      probs = matrix(0, nrow=S, ncol=G_max)
      means = array(rep(0, S*max(J)*G_max), dim=c(S,max(J),G_max))
      covs = array(rep(0, S*max(J)*max(J)*G_max), dim=c(S,max(J),max(J),G_max))
      for (s in 1:S) {
        for (j in 1:J[s]) {
          logblhaz[, s, j] =
            smpl.df[, paste0('logblhaz[', s, ',', j, ']')] # skip the normal approx step
          # rnorm(n = nrow(smpl.df),
          #       mean = smpl.df[, paste0('logblhaz_mean[', s, ']')],
          #       sd = smpl.df[, paste0('logblhaz_sd[', s, ']')])
        }
        bhm_approx = Mclust(logblhaz[,s,1:J[s]], G=1:G_max)
        G[s] = bhm_approx$G
        probs[s, 1:G[s]] = bhm_approx$parameters$pro
        means[s, 1:J[s], 1:G[s]] = bhm_approx$parameters$mean
        covs[s, 1:J[s], 1:J[s], 1:G[s]] = bhm_approx$parameters$variance$sigma
      }
      save(list=c('G', 'probs', 'means', 'covs'),
           file=paste0("data/bhm_sanofi_", ifelse(gene==1, 'wt', 'mutant'), ifelse(run.ps, '_ps', ''), ".Rdata"))
    } else {
      load(paste0("data/bhm_sanofi_", ifelse(gene==1, 'wt', 'mutant'), ifelse(run.ps, '_ps', ''), ".Rdata"))
      G_max = 5
      
      # 3. RMAP portion
      model_rmapp = cmdstan_model('rmapp_pwe.stan')
      smpl_rmap <- model_rmapp$sample(
        data = list(
          'S'=S, 'J'=J,
          'G'=G, 'G_max'=G_max, 'probs'=probs, 'means'=means, 'covs'=covs,
          'w'=w,
          'norm_vague_mean'=0, 'norm_vague_sd'=10,
          'n' = n, 'Ji' = Ji, 'event' = event.long, 'logrisktime' = logrisktime.long, 'x' = x,
          'beta_mean' = 0, 'beta_sd' = 5
          # 'beta_mean_mean' = 0, 'beta_mean_sd' = 10,
          # 'beta_sd_scale' = 1, 'beta_sd_df' = 25
        )
        , iter_warmup = 2000, iter_sampling = nsample+2000, chains = 1, parallel_chains = 1
        , init=function() list(#blhaz=ifelse(is.na(h), 0, h),
          beta=rep(0, S))
      )
      
      smpl <- posterior::as_draws_df(smpl_rmap$draws()) %>% as.data.frame() %>%
        select(starts_with("blhaz") | starts_with("beta"))
      # colMeans(smpl)%>% head(30) %>% exp()
      # hist(exp(smpl$`logblhaz[1,1]`))
    }
  } else if (method=='cwpp') {
    #------------------#
    # Case Weighted PP (Kwiatkowski et al)
    #------------------
    histdata=a
    curdata=b
    histcur=c
    res_all_0 = data.frame()
    for (s in 1:S) {
      sim.data = curdata[,c('os', 'death', 'trt', 's')] %>%
        rename(t=os,
               nu=death,
               treat=trt,
               strata=s) %>%
        filter(strata==s) %>%
        mutate(external=0, nu=as.numeric(nu)) %>%
        rbind(
          histdata[,c('os', 'death', 'trt', 's')] %>%
            rename(t=os,
                   nu=death,
                   treat=trt,
                   strata=s) %>%
            filter(strata==s) %>%
            mutate(external=1, nu=as.numeric(nu)))
      if (sum(sim.data$external==1)==0 ||
        !any(sim.data$nu[sim.data$external==1]==0)) next # skip if historical has no censoring (the cwpp code can't handle that)
      source('/proj/ibrahimlab/leap/sims_pwe/case_weighted_pp.R')
      res_all_0 = res_all_0 %>% rbind(res_s%>%t%>%as.data.frame%>%mutate(strata=s))
    }
    a=histdata
    b=curdata
    c=histcur
    if(any(grepl("package:MASS", search())))
      detach("package:MASS", unload=TRUE) # remove MASS library or dplyr won't work..
    
    res_all_0 = res_all_0 %>%
      mutate(variable_beta = paste0('beta[', strata, ']'),
             variable_HR = paste0('expbeta[', strata, ']'))
    res_all_beta = res_all_0 %>%
      select(variable_beta, a0_logHR, SE_a0, lower_logHR_a0, upper_logHR_a0) %>%
      rename(variable = variable_beta,
             mean = a0_logHR,
             sd = SE_a0,
             `2.5%` = lower_logHR_a0,
             `97.5%` = upper_logHR_a0)
    res_all_HR = res_all_0 %>%
      select(variable_HR, HR_a0, HR_SE_a0, lower_a0, upper_a0) %>%
      rename(variable = variable_HR,
             mean = HR_a0,
             sd = HR_SE_a0,
             `2.5%` = lower_a0,
             `97.5%` = upper_a0)
    res[[method]] = res_all_beta %>% rbind(res_all_HR) %>%
      mutate(median=NA,
             ess_bulk=NA,
             method='cwpp')
  }
  
  if (runBHMonly) break
  if (pooled) method='bhm_pooled'
  
  if (method=='noninf') {
    beta_mean = smpl[,sapply(1:S, function(s) paste0('beta[', s, ']'))] %>% rowMeans
    beta_sd = smpl[,sapply(1:S, function(s) paste0('beta[', s, ']'))] %>% apply(1, sd)
    smpl = cbind(beta_mean, beta_sd, smpl)
  }
  
  if (method!='cwpp') {
    smpl.mcmc <- as_draws_matrix(smpl)[-1,]
    # smpl.mcmc %>% mcmc_trace('beta[1]')
    # smpl.mcmc %>% mcmc_hist('beta[1]')
    res[[method]] = 
      smpl.mcmc %>% summarize_draws(
        'mean', 'median', 'sd', ~quantile(.x, probs = c(0.025, 0.975))
        , 'ess_bulk'
      ) %>%
      mutate(across(where(is.numeric), round, 3),
             method=method)
  }
  

  smpls[[method]] = smpl.mcmc
}

## Fit Cox
if (gene==1) {
  coxfit   <- lapply(1:S, function(x) coxph(Surv(os, death) ~ trt, data = b %>% filter(s == x)))
} else { # last strata has 1 obs so do not run cox
  coxfit   <- lapply(1:(S-1), function(x) coxph(Surv(os, death) ~ trt, data = b %>% filter(s == x)))
}
coxfit[[length(coxfit)+1]] = coxph(Surv(os, death) ~ trt + strata(s), data = b)
coxmle   <- lapply(coxfit, function(f) summary(f)$coefficients[c(1, 3)])
cox.summ <- t(sapply(coxmle, function(f) c(f, f[1] + c(-1, 1) * qnorm(0.975) * f[2])))
# get exp(beta) version
coxmle2   <- lapply(coxfit, function(f) c(summary(f)$conf.int[1], # exp(beta) mean
                                          summary(f)$conf.int[1]*summary(f)$coefficients[3], # SE
                                          summary(f)$conf.int[3:4])) # CI
cox.summ = rbind(cox.summ, do.call(rbind, coxmle2))
if (gene==2) cox.summ = rbind(cox.summ[1:7,],
                              rep(NA, 4),
                              cox.summ[8,],
                              cox.summ[9:15,],
                              rep(NA, 4),
                              cox.summ[16,])
cox.summ <- data.frame(variable = c(sapply(1:S, function(i) paste0('beta[',i,']')), 'beta_mean',
                                    sapply(1:S, function(i) paste0('expbeta[',i,']')), 'expbeta_mean'),
                       cox.summ,
                       ess_bulk=NA,
                       method = 'cox')
names(cox.summ) <- c('variable', 'mean', 'sd', '2.5%', '97.5%', 'ess_bulk', 'method')
res[['cox']] = cox.summ %>% mutate(median=NA)

res = do.call('rbind', res)
saveRDS(res, paste0('data_analysis_res_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.rds'))
save(smpls, file = paste0('data_analysis_smpls_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.RData'))

#--------------------#

res = list()#res = readRDS(paste0('data_analysis_res_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.rds'))
load(paste0('data_analysis_smpls_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.RData'))

## Tranform beta (log HR) to HR
for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
  smpls[[method]] = smpls[[method]] %>%
    as.data.frame %>%
    mutate(across(starts_with("beta["), exp, .names = "exp{.col}"),
           expbeta_mean = exp(beta_mean)) %>%
    as_draws_matrix
  
  res[[method]] = smpls[[method]] %>%
    summarize_draws(
      'mean', 'median', 'sd', ~quantile(.x, probs = c(0.025, 0.975))
      , 'ess_bulk'
    ) %>%
    mutate(across(where(is.numeric), round, 3),
           method=method)
}
res = do.call('rbind', res)
res = res %>% rbind(
  readRDS(paste0('data_analysis_res_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.rds')) %>%
    filter(method%in%c('cox','cwpp'))
)

# ---------------------
## sLEAP Probabilities
# ---------------------

png(paste0('data_analysis_prob_', ifelse(gene==1, 'WT', 'MUTANT'), ifelse(run.ps, '_ps', ''), '.png'),
    width = 5000, height = 2000, res=600)
ggplot(res %>% filter(grepl('^gamma', variable)) %>% 
         mutate(s = case_when(variable == 'gamma[1, 1]' ~  1,
                              variable == 'gamma[2, 1]' ~  2,
                              variable == 'gamma[3, 1]' ~  3,
                              variable == 'gamma[4, 1]' ~  4,
                              variable == 'gamma[5, 1]' ~  5,
                              variable == 'gamma[6, 1]' ~  6,
                              variable == 'gamma[7, 1]' ~  7,
                              variable == 'gamma[8, 1]' ~  8)) %>%
         filter(grepl(', 1]', variable)),
       aes(s, mean)) +
  geom_line(linewidth=0.8, alpha=.6) +
  theme_minimal() +
  geom_errorbar(aes(ymin=`2.5%`,
                    ymax=`97.5%`), width=.1) +
  scale_x_continuous(breaks = 1:8) +
  labs(y = 'mean post. gamma, 95% CI', x = 'Stratum')
dev.off()

# --------------
## Get DIC
# --------------

# Calculate PWE Log-Likelihood using long data
get_loglik <- function(d, beta=NULL, lambda=NULL) {
  if (is.null(beta)) { # Under observed data.. this might be wrong
    hazard = d$death
  } else {
    hazard = lambda[cbind(d$s, d$intervalind)] * exp(d$trt * beta[d$s])
  }
  loglik = d$death * log(hazard) + hazard * d$risktime
  loglik[is.nan(loglik)] = 0
  sum(loglik)
}

d = c.long %>% filter(histdata==0)

dic = list()
for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
  smpl = smpls[[method]]
  
  # Get posterior mean of log likelihoods
  post.logliks = rep(NA, nrow(smpl))
  for (m in 1:nrow(smpl)) {
    beta = c(smpl[m, sapply(1:S, function(s) paste0('beta[', s, ']'))])
    lambda = matrix(nrow=S, ncol=max(J))
    for (s in 1:S) {
      for (j in 1:J[s]) {
        if (method=='sleap') {
          lambda[s, j] = smpl[m, paste0('blhaz[', s, ', 1, ', j, ']')]
        } else {
          lambda[s, j] = smpl[m, paste0('blhaz[', s, ', ', j, ']')]
        }
      }
    }
    post.logliks[m] = get_loglik(d, beta, lambda)
  }
  postmean.dev =  -2 * mean(post.logliks)
  
  # Get log likelihoods under posterior means
  postmeans = smpl %>% colMeans
  postmean.beta = postmeans[sapply(1:S, function(s) paste0('beta[', s, ']'))]
  postmean.lambda = matrix(nrow=S, ncol=max(J))
  for (s in 1:S) {
    for (j in 1:J[s]) {
      if (method=='sleap') {
        postmean.lambda[s, j] = postmeans[paste0('blhaz[', s, ', 1, ', j, ']')]
      } else {
        postmean.lambda[s, j] = postmeans[paste0('blhaz[', s, ', ', j, ']')]
      }
    }
  }
  dev.postmean = -2 * get_loglik(d, postmean.beta, postmean.lambda)
  
  pd = postmean.dev - dev.postmean
  dic[[method]] = pd + postmean.dev
}
dic[['cox']] <- dic[['cwpp']] <- NA

# --------------
## Print latex
# --------------

res %>%
  filter(variable=='expbeta_mean') %>%
  mutate(ci = paste0(formatC(`2.5%`, format = 'f', digits = 3), ', ', formatC(`97.5%`, format = 'f', digits = 3)),
         method = factor(method, levels=c('sleap', 'pp', 'cwpp', 'bhm', 'noninf', 'cox')),
         dic = as.numeric(dic[method]),
         dic.char = ifelse(is.na(dic), '', formatC(dic, format = 'f', digits = 2)),
         mean = formatC(mean, format = 'f', digits = 3),
         sd = formatC(sd, format = 'f', digits = 3)) %>%
  select(method, variable, mean, sd, ci, dic.char) %>%
  tabular(Factor(method, 'Prior', levelnames = c('sLEAP + BHM', 'NPP + BHM', 'CWPP', 'BHM', 'NonInf', 'Cox')) ~
            Heading() * identity * (
              Heading('DIC') * dic.char +
                Heading('Treatment') * (
                  Heading('Mean') * mean + Heading('SD') * sd + Heading('95\\% CI') * ci
                )
            ), data=.) %>%
  toLatex %>%
  print

# --------------------------------------------------------
## Print latex for stratum table and fix strata ordering
# --------------------------------------------------------

## Create strata key-table
strata.table = c %>%
  group_by(ecog2, bevacizumab, oxaliplatin) %>%
  summarize(n0=sum(histdata==1), n=sum(histdata==0)) %>%
  arrange(desc(n0)) %>%
  ungroup() %>%
  mutate(s=row_number())

res.betas = res %>%
  filter(grepl('expbeta\\[', variable),
         method=='sleap') %>%
  cbind(J)

res.betas2 = res.betas[match(strata.table$n, n),] %>% # CHANGE THE ORDER OF STRATA TO MATCH THE STRATA TABLE
  cbind(strata.table) %>%
  mutate(strata = paste0('(', ecog2, ', ', bevacizumab, ', ', oxaliplatin, ')'),
         ci = paste0(formatC(`2.5%`, format = 'f', digits = 3), ', ', formatC(`97.5%`, format = 'f', digits = 3)),
         mean = formatC(mean, format = 'f', digits = 3),
         sd = formatC(sd, format = 'f', digits = 3)) %>%
  select(strata, n0, n, J, mean, sd, ci)

res.betas2 %>%
  mutate(strata = factor(strata, levels=res.betas2$strata)) %>%
  tabular(Factor(strata, 'Strata') ~
            Heading() * identity * (
              Heading('\\|J\\|') * J +
                Heading('n_0') * n0 +
                Heading('n') * n +
                Heading('Treatment') * (
                  Heading('Mean') * mean + Heading('SD') * sd + Heading('95\\% CI') * ci
                )
            ), data=.) %>%
  toLatex %>%
  print
