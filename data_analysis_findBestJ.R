library(dplyr)
library(tidyr)
library(survival)
library(nimble)
library(posterior)
library(bayesplot)
library(cmdstanr)
library(mclust)
library(ggplot2)

gene = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
nsample = 10000
runBHMonly = F
pooled = F

source('models.R')

a = read.csv('data/surv_sanofi.csv') %>% rename_all(tolower) %>% rename(subjid=rsubjid) %>%
  mutate(histdata=1, kras=NA, race=case_when(race=='CAUCASIAN/WHITE' ~ 'White or Caucasian',
                                             race=='BLACK' ~ 'Black',
                                             race=='OTHER' ~ 'Other'))
b = read.csv('data/surv_amgen.csv') %>% rename_all(tolower) %>% mutate(histdata=0) %>%
  filter(kras==ifelse(gene==1, 'Wild-type', 'Mutant')) %>%
  # filter(pfs!=0) %>% # only use for PFS
  mutate(kras=ifelse(kras=='', 'NA', kras))
c = bind_rows(a, b) %>%
  mutate(across(all_of(c('ecog', 'ecog2', 'bevacizumab', 'oxaliplatin', 'histdata', 'kras', 'race', 'sex')), as.factor),
         # trt = factor(trt, level=c('1', '0')), # comment out if not doing descriptive stats.
         s = as.numeric(interaction(ecog2, bevacizumab, oxaliplatin)))
# interaction(c$ecog2, c$bevacizumab, c$oxaliplatin)['Levels'] # see strata numbers
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
n0 = a %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)
n = b %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)

# -----------------------------
# FIND BEST J FOR EACH STRATUM
# -----------------------------

a.orig = a
b.orig = b
c.orig = c
S.orig = S
n0.orig = n0
n.orig = n
J.max = ifelse(ceiling(n.orig/5)>10, 10, ceiling(n.orig/5))

dics = matrix(nrow=S.orig, ncol=max(J.max))

for (s.orig in 1:S.orig) {
  
  if (J.max[s.orig]==1) next
  
  a = a.orig %>% filter(s==s.orig) %>% mutate(s=1)
  b = b.orig %>% filter(s==s.orig) %>% mutate(s=1)
  c = c.orig %>% filter(s==s.orig) %>% mutate(s=1)
  
  S=max(c$s)
  K=rep(2, S)
  
  n0 = a %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)
  n = b %>% group_by(s) %>% summarize(n=n()) %>% left_join(tibble(s=1:max(c$s)),., by='s') %>% mutate(n=ifelse(is.na(n), 0, n)) %>% pull(n)
  
  for (J in 1:J.max[s.orig]) {
    
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
    
    #------------------#
    # NON-INFORMATIVE
    #------------------
    niminit  <- list(
      beta = setup$inits$beta
      , blhaz = setup$inits$blhaz[,1,] %>% matrix(nrow=1)
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
    
    smpl <- as_draws_matrix(smpl)[-1,]
    
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
    
    # Get posterior mean of log likelihoods
    post.logliks = rep(NA, nrow(smpl))
    for (m in 1:nrow(smpl)) {
      beta = c(smpl[m, sapply(1:S, function(s) paste0('beta[', s, ']'))])
      lambda = matrix(nrow=S, ncol=max(J))
      for (s in 1:S) {
        for (j in 1:J[s]) {
          lambda[s, j] = smpl[m, paste0('blhaz[', s, ', ', j, ']')]
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
        postmean.lambda[s, j] = postmeans[paste0('blhaz[', s, ', ', j, ']')]
      }
    }
    dev.postmean = -2 * get_loglik(d, postmean.beta, postmean.lambda)
    
    pd = postmean.dev - dev.postmean
    dics[s.orig, J] = pd + postmean.dev
  }
}

dics
saveRDS(dics, paste0('data_analysis_findBestJ_', ifelse(gene==1, 'WT', 'MUTANT'), '.rds'))
