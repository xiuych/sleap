library(dplyr)
library(optimx)
library(numDeriv)
library(survival)
library(ggplot2)
library(ggfortify)
remove(list = ls())

source("helper_fun.R")
load("data/mle.Rdata")

grid = expand.grid(p = c(0.0, 0.5, 1.0),
                   n = c(25, 50, 100),
                   q = c(-0.4, -0.2, 0, 0.2, 0.4))

id=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#-----------------------------
# Generate historical data
#-----------------------------

histdata.all = data.frame()
for (s in 1:S) {
  
  load("data/mle.Rdata")
  n                 <- grid$n[id]     ## number of subjects in current data
  n0                <- grid$n[id]     ## number of subjects in historical data
  prob.trt          <- trt.prop## probability of receiving treatment (current data)
  prob.exch         <- grid$p[id]     ## probability of being exchangeable (historical data)
  q                 <- grid$q[id]     ## unexchangeability parameter: Y0 | K = k ~ Weibull( shape, scale * exp(q * k ) )
  # q positive ->scale param increase ->survival decrease in nonexch group
  shape             <- blshape[s]   ## shape parameter for Weibull distribution
  scale             <- scale[s]     ## scale parameter for Weibull distribution
  beta              <- trteff[s]    ## regression coefficient: Y | A = a ~ Weibull(shape, scale * exp(beta * a))
  # beta positive -> scale param increase -> surv decrease in trt group i.e. control group healthier
  prop.event.ctrl   <- event.rate.ctrl[s]    ## proportion of observed event times, ctrl group
  prop.event.trt    <- event.rate.trt[s]     ## proportion of observed event times, trt group
  hist.thresh       <- 0.01    ## need: |mle - truth| < hist.thresh using [log(shape), log(scale), q]
  hist.max          <- 100000  ## max number of iterations to generate historical data
  
  
  ## GENERATE HISTORICAL DATA VIA LOOP
  if ( prob.exch == 1 )
    truth <- log( c(shape, scale) )
  if ( prob.exch == 0 )
    truth <- log( c(shape, exp(q) * scale) )
  if ( !(prob.exch %in% c(0, 1) ) )
    truth <- c(log(shape), log(scale), q)
  criterion  <- 100
  best.seed  <- 0
  best.mle   <- NULL
  
  pb <- txtProgressBar(
    min = 0, max = hist.max, style = 3, width = 50, char = "="
  )
  
  for ( i in 1:hist.max ) {
    setTxtProgressBar(pb, i)
    set.seed(i)
    ## Generate historical data using prob.trt = prob.exch and beta = q
    histdata           <- gendata(n0, prob.trt = 1-prob.exch, shape, scale, beta = q, prop.event.trt = prop.event.ctrl, prop.event.ctrl = prop.event.ctrl)
    names(histdata)[3] <- 'k'
    ## Obtain correct function to optimize; if prob.exch == 0 or 1, there are only 2 parameters
    if ( prob.exch == 0 | prob.exch == 1 ) {
      fn <- function(parms) negloglik(parms, histdata$y, histdata$event, histdata$k, ctrlonly = TRUE)
      gr <- function(parms) gnegloglik(parms, histdata$y, histdata$event, histdata$k, ctrlonly = TRUE)
    } else {
      fn <- function(parms) negloglik(parms, histdata$y, histdata$event, histdata$k, ctrlonly = FALSE)
      gr <- function(parms) gnegloglik(parms, histdata$y, histdata$event, histdata$k, ctrlonly = FALSE)
    }
    ## Compute MLE for historical data using non-exchangeability indicator as treatment indicator.
    ## The estimated treatment effect should be ~= q ideally
    opt <- optim(truth, fn = fn, gr = gr, method = 'BFGS')
    mle        <- opt$par
    ## Check if the MLE for the historical data set is close to (log(shape), log(scale), q)
    prev       <- criterion
    criterion  <- max(abs(truth - mle))
    if ( criterion < prev ) {
      best.seed <- i
      best.mle  <- mle
    }
    if ( criterion < hist.thresh ) {
      cat('\n')
      print(paste0('Converged at seed = ', best.seed))
      break
    }
  }
  ## Provide warning if there was no convergence; report best possible seed
  ## and re-generate the historical data set using that seed.
  if ( i == hist.max ) {
    cat('\n')
    warning(paste0('Did not converge. Best seed is ', best.seed))
    set.seed(best.seed)
    histdata <- gendata(n0, prob.trt = prob.exch, shape, scale, beta = q, prop.event.trt = prop.event.ctrl, prop.event.ctrl = prop.event.ctrl)
    names(histdata)[3] <- 'k'
  }
  histdata$a    <- 0
  histdata$exch <- with(histdata, 1 - k)
  
  histdata.all = histdata %>% mutate(strata=s) %>% rbind(histdata.all)
}

saveRDS(histdata.all, file=paste0('data/historical/histdata_', id, '.rds'))

#--------------------------
# Generate current data
#--------------------------
# 
# for (n_gen in c(25, 50, 100)) {
#   for (i in 1:10000) {
#     if (!is.na(id)) break
# 
#     data.all = data.frame()
#     for (s in 1:S) {
# 
#       load("data/mle.Rdata")
#       n0 <- n <- n_gen
#       prob.trt          <- trt.prop## probability of receiving treatment (current data)
#       shape             <- blshape[s]   ## shape parameter for Weibull distribution
#       scale             <- scale[s]     ## scale parameter for Weibull distribution
#       beta              <- trteff[s]    ## regression coefficient: Y | A = a ~ Weibull(shape, scale * exp(beta * a))
#       prop.event.ctrl   <- event.rate.ctrl[s]    ## proportion of observed event times, ctrl group
#       prop.event.trt    <- event.rate.trt[s]     ## proportion of observed event times, trt group
#       hist.thresh       <- 0.01    ## need: |mle - truth| < hist.thresh using [log(shape), log(scale), q]
#       hist.max          <- 100000  ## max number of iterations to generate historical data
# 
#       ## GENERATE CURRENT DATA
#       data      <- gendata(n, prob.trt, shape, scale, beta, prop.event.trt = prop.event.trt, prop.event.ctrl = prop.event.ctrl)
#       data$exch <- 1
# 
#       # attributes(histdata)$params[c('cens.rate.trt', 'cens.rate.ctrl')]
#       # attributes(data)$params[c('cens.rate.trt', 'cens.rate.ctrl')]
#       #
#       # ## Fit KM curve to see similarity
#       # merge     <- histdata %>% mutate(data = 'hist') %>% merge(y = data %>% mutate(data = 'current'), all = TRUE) %>%
#       #   mutate(trt = ifelse(a == 1, 'trt', 'ctrl'), ex = ifelse(exch == 1, 'exch', 'unexch'))
#       # model_fit <- survfit(Surv(y, event) ~  data + exch, data = merge %>% filter(a==0))
#       # autoplot(model_fit, se = FALSE)
#       #
#       # cbind(
#       #   'hist.exch'     = summary(histdata %>% filter(exch == 1) %>% select(t))
#       #   , 'cur.ctrl'    = summary(data %>% filter(a == 0) %>% select(t))
#       #   , 'hist.unexch' = summary(histdata %>% filter(exch == 0) %>% select(t))
#       # )
# 
#       data.all = data %>% mutate(strata=s) %>% rbind(data.all)
#     }
# 
#     saveRDS(data.all, file=paste0('data/current/curdata_', n, '_', i, '.rds'))
#   }
# }
