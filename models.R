sleap_wide_data <- function (
    y, event, breaks
) {
  n <- length(y)
  J <- length(breaks) - 1
  logrisktime <- matrix(-100, nrow = n, ncol = J)
  eventind    <- matrix(-100, nrow = n, ncol = J)
  Ji          <- rep(-1, n)
  colnames(logrisktime) <- paste0('logrt_', 1:J)
  colnames(eventind) <- paste0('event_', 1:J)
  for ( j in 1:J ) {
    indx_action   <- ( y > breaks[j] ) & ( y <= breaks[j+1] )
    indx_noaction <- ( y > breaks[j+1] )
    eventind[indx_action, j]    <- event[indx_action]
    eventind[indx_noaction, j]  <- 0
    logrisktime[indx_action, j]   <- log( y[indx_action] - breaks[j])
    logrisktime[indx_noaction, j] <- log( breaks[j+1] - breaks[j] )
    Ji[indx_action]  <- j
  }
  data.frame(eventind, logrisktime, Ji)
}

pwe.mle <- function(
    widedata, xnames = NULL
) {
  J <- max(widedata$Ji)
  longdat <- widedata %>% 
    pivot_longer(
      -c("Ji", any_of(xnames)), names_to = c(".value", "interval"), names_pattern = "(.+)_(.+)"
    ) %>%
    filter(event >= 0)
  
  if ( J > 1 ) {
    pwe.formula <- as.formula( paste0( c( 'event ~ 0 + offset(logrt) + factor(interval)', xnames ), collapse = ' + ' ) )
  } else {
    pwe.formula <- as.formula( paste0( c( 'event ~ 1 + offset(logrt)', xnames ), collapse = ' + ' ) )
  }
  glm(pwe.formula, 'poisson', longdat)
}

sleap_setup <- function(
    formula, data, histdata, stratavar, nclasses
    , nintervals = NULL, pooled_intervals = FALSE
    , bhm = TRUE
    , beta_mean = 0, beta_sd = 5
    , beta_mean_mean = 0, beta_mean_sd = 10
    , beta_sd_mu = 0, beta_sd_sigma = 0.5
    , blhaz_shape = 0.1, blhaz_rate = 0.1
    , conc = 0.95
) {
  
  ## Create pooled data set
  data$hist     <- 0
  histdata$hist <- 1
  pooled        <- rbind(data, histdata)
  
  ## Get names based on formula
  varnames  <- all.vars(formula)
  timename  <- varnames[1]
  eventname <- varnames[2]
  xnames    <- varnames[-(1:2)]
  
  ## Drop any observations whose strata ID is not in current data
  strat.cur.unique    <- sort( unique(data[[stratavar]]) )
  nstrata             <- length(strat.cur.unique)
  histdata            <- histdata %>% filter(get(stratavar) %in% strat.cur.unique)
  pooled[[stratavar]] <- factor(pooled[[stratavar]])
  pooled$id           <- 1:nrow(pooled)
  
  Jstar <- max(nintervals)
  Kstar <- max(nclasses)
  
  ## Loop through strata; create stratum-specific breaks; compute time at risk
  data.list   <- vector('list', nstrata)
  breaks.list <- data.list
  names(data.list) <- names(breaks.list) <- strat.cur.unique
  
  data.list        <- vector('list', nstrata)
  names(data.list) <- strat.cur.unique
  breaks.list <- fit.cur <- fit.hist <- data.list
  beta.init   <- rep(NA, nstrata)
  gamma.init  <- array(0, dim = c(nstrata, Kstar))
  blhaz.init  <- array(-100, dim = c(nstrata, Kstar, Jstar))
  c0.init     <- sample(1:2, nrow(histdata), replace = TRUE)
  ## Loop through strata; create stratum-specific breaks; compute time at risk
  for ( s in 1:nstrata ) {
    gamma.init[s, 1:nclasses[s]] <- 1 / nclasses[s]
    ## Create intervals based on observed event times in stratum
    data.list[[s]] <- pooled %>% filter(get(stratavar) == strat.cur.unique[s])
    if ( pooled_intervals ) {
      eventtimes  <- data.list[[s]] %>% filter(get(eventname) == 1) %>% ungroup() %>% select(all_of(timename)) %>% unlist
    } else {
      eventtimes  <- data.list[[s]] %>% filter(get(eventname) == 1, hist == 0) %>% ungroup() %>% select(all_of(timename)) %>% unlist
    }
    if ( nintervals[s] > 1 ) {
      breaks.list[[s]] <- quantile(eventtimes, probs = (1:(nintervals[s]-1)) / (nintervals[s]))
      breaks.list[[s]] <- c(0, breaks.list[[s]], Inf)
    } else {
      breaks.list[[s]] <- c(0, Inf)
    }
    remove(list = 'eventtimes')
    
    ## Obtain wide data for Nimble
    temp <- sleap_wide_data(
      y = data.list[[s]][[timename]]
      , event = data.list[[s]][[eventname]]
      , breaks = breaks.list[[s]]
    )
    data.list[[s]]         <- cbind(temp, data.list[[s]] %>% select(xnames, hist))
    data.list[[s]]$stratum <- s
    
    ## Fit MLE to current / historical to obtain starting values
    fit.cur[[s]]       <- pwe.mle(data.list[[s]] %>% filter(hist == 0), xnames = xnames)
    if (data.list[[s]] %>% filter(hist == 1) %>% nrow() != 0)
      fit.hist[[s]]      <- pwe.mle(data.list[[s]] %>% filter(hist == 1), xnames = NULL)
    beta.init[s]  <- coef(fit.cur[[s]])[xnames]
    blhaz.init[s, 1, 1:nintervals[s]] <- exp( coef(fit.cur[[s]])[1:nintervals[s]] )
    if (data.list[[s]] %>% filter(hist == 1) %>% nrow() != 0)
      blhaz.init[s, 2, 1:nintervals[s]] <- exp( coef(fit.hist[[s]])[1:nintervals[s]] )
  }
  ## Stack all data--fill in NA where necessary
  data.stacked <- plyr::rbind.fill(data.list)
  ## Replace NA with nonmissing
  data.stacked <- data.stacked %>% replace(is.na(.), -100)
  
  
  res <- list(
    const = list(
      n = data.stacked %>% filter(hist == 0) %>% nrow
      , n0 = data.stacked %>% filter(hist == 1) %>% nrow
      , J = nintervals, S = nstrata, K = nclasses
      , Ji        = data.stacked %>% filter(hist == 0) %>% select('Ji') %>% unlist
      , Ji0          = data.stacked %>% filter(hist == 1) %>% select('Ji') %>% unlist
      , strat     = as.numeric( data.stacked %>% filter(hist == 0) %>% select(stratum) %>% unlist )
      , strat0       = as.numeric( data.stacked %>% filter(hist == 1) %>% select(stratum) %>% unlist )
    )
    , data = list(
      logrisktime = data.stacked %>% filter(hist == 0) %>% select(starts_with('logrt_')) %>% as.matrix
      , event     = data.stacked %>% filter(hist == 0) %>% select(starts_with('event_')) %>% as.matrix
      , logrisktime0 = data.stacked %>% filter(hist == 1) %>% select(starts_with('logrt_')) %>% as.matrix
      , event0       = data.stacked %>% filter(hist == 1) %>% select(starts_with('event_')) %>% as.matrix
      , x              = data.stacked %>% filter(hist == 0) %>% select(all_of(xnames)) %>% unlist
      , blhaz_rate     = blhaz_rate
      , blhaz_shape    = blhaz_shape
      , conc           = rep(conc, max(nclasses))
      , beta_sd_mu     = beta_sd_mu
      , beta_sd_sigma  = beta_sd_sigma
    )
    , inits = list(
      beta = beta.init, blhaz = blhaz.init, gamma = gamma.init, c0 = c0.init
    )
    , rawdata = data.stacked
    , mle.cur = fit.cur
    , mle.hist = fit.hist
  )
  if (bhm) {
    res$inits$beta_mean = 0
    res$inits$beta_sd = 1
    res$data$beta_mean_mean = beta_mean_mean
    res$data$beta_mean_sd   = beta_mean_sd
    res$data$beta_sd_mu     = beta_sd_mu
    res$data$beta_sd_sigma  = beta_sd_sigma
    
  } else {
    res$data$beta_mean = beta_mean
    res$data$beta_sd   = beta_sd
  }
  res
}

sleap <- nimbleCode({
  beta_mean ~ dnorm(beta_mean_mean, sd = beta_mean_sd)
  beta_sd   ~ T( dnorm(beta_sd_mu, sd = beta_sd_sigma), 0, Inf )
  for ( s in 1:S ) {
    ## Normal priors on treatment effects
    beta[s] ~ dnorm(mean = beta_mean, sd = beta_sd)
    ## Dirichlet priors on exchangeability parameters
    gamma[s,1:K[s]] ~ ddirch(conc[1:K[s]])
    ## Gamma priors on baseline hazards
    for ( k in 1:K[s] ) {
      for ( j in 1:J[s] ) {
        blhaz[s, k, j] ~ dgamma(shape = blhaz_shape, rate = blhaz_rate)
      }
    }
  }
  ## Likelihood of current data
  for ( i in 1:n ) {
    eta[i] <- x[i] * beta[strat[i]]
    for ( j in 1:Ji[i] ) {
      event[i, j] ~ dpois( blhaz[strat[i], 1, j] * exp( eta[i] + logrisktime[i,j] ) )
    }
  }
  ## Complete data likelihood of historical data
  for ( i in 1:n0 ) {
    c0[i] ~ dcat(prob = gamma[strat0[i], 1:K[strat0[i]]])
    for ( j in 1:Ji0[i] ) {
      event0[i, j] ~ dpois( blhaz[strat0[i], c0[i], j] * exp( logrisktime0[i,j] ) )
    }
  }
})

sleap2 <- nimbleCode({
  for ( s in 1:S ) {
    ## Normal priors on treatment effects
    beta[s] ~ dnorm(mean = beta_mean, sd = beta_sd)
    ## Dirichlet priors on exchangeability parameters
    gamma[s,1:K[s]] ~ ddirch(conc[1:K[s]])
    ## Gamma priors on baseline hazards
    for ( k in 1:K[s] ) {
      for ( j in 1:J[s] ) {
        blhaz[s, k, j] ~ dgamma(shape = blhaz_shape, rate = blhaz_rate)
      }
    }
  }
  ## Likelihood of current data
  for ( i in 1:n ) {
    eta[i] <- x[i] * beta[strat[i]]
    for ( j in 1:Ji[i] ) {
      event[i, j] ~ dpois( blhaz[strat[i], 1, j] * exp( eta[i] + logrisktime[i,j] ) )
    }
  }
  ## Complete data likelihood of historical data
  for ( i in 1:n0 ) {
    c0[i] ~ dcat(prob = gamma[strat0[i], 1:K[strat0[i]]])
    for ( j in 1:Ji0[i] ) {
      event0[i, j] ~ dpois( blhaz[strat0[i], c0[i], j] * exp( logrisktime0[i,j] ) )
    }
  }
})

## Power Prior code (stratum-specific)
pp <- nimbleCode({
  ## Hyperprior for treatment effects
  beta_mean ~ dnorm(beta_mean_mean, sd = beta_mean_sd)
  beta_sd   ~T ( dnorm(beta_sd_mu, sd = beta_sd_sigma), 0, Inf )
  
  ## LOOP THROUGH NUMBER OF STRATA
  for ( s in 1:S ) {
    a0[s] ~ dunif(0, 1)
    beta[s] ~ dnorm(beta_mean, sd = beta_sd)
    
    ## LIKELIHOOD OF CURRENT DATA
    for ( i in 1:n[s] ) {
      linpred[s,i] <- x[s,i] * beta[s]
      for ( j in 1:Ji[s,i] )
        event[s,i,j] ~ dpois( blhaz[s, j] * exp( linpred[s,i] + logrisktime[s,i,j] ) )
    }
    
    ## Conjugate prior using historical data
    for (j in 1:J[s]) {
      blhaz[s,j] ~ dgamma(blhaz_shape + a0[s] * sum(event0[s,,j]),
                          blhaz_rate + a0[s] * sum(exp(logrisktime0[s,,j])))
    }
  }
})

## Power Prior code (stratum-specific) (no bhm)
pp2 <- nimbleCode({
  ## LOOP THROUGH NUMBER OF STRATA
  for ( s in 1:S ) {
    a0[s] ~ dunif(0, 1)
    beta[s] ~ dnorm(beta_mean, sd = beta_sd)
    
    ## LIKELIHOOD OF CURRENT DATA
    for ( i in 1:n[s] ) {
      linpred[s,i] <- x[s,i] * beta[s]
      for ( j in 1:Ji[s,i] )
        event[s,i,j] ~ dpois( blhaz[s, j] * exp( linpred[s,i] + logrisktime[s,i,j] ) )
    }
    
    ## Conjugate prior using historical data
    for (j in 1:J[s]) {
      blhaz[s,j] ~ dgamma(blhaz_shape + a0[s] * sum(event0[s,,j]),
                          blhaz_rate + a0[s] * sum(exp(logrisktime0[s,,j])))
    }
  }
})

## BHM prior code
bhm <- nimbleCode({
  ## Hyperprior for exchangeable treatment effects
  beta_mean ~ dnorm(beta_mean_mean, sd = beta_mean_sd)
  beta_sd   ~ T( dnorm(beta_sd_mu, sd = beta_sd_sigma), 0, Inf )
  
  ## LOOP THROUGH NUMBER OF STRATA
  for ( s in 1:S ) {
    ## PRIORS FOR EXCHANGEABLE PARAMETERS
    ## Prior for exchangeable treatment effects
    beta[s] ~ dnorm(beta_mean, sd = beta_sd)
    for ( j in 1:J[s] )
      blhaz[s,j] ~ dgamma(blhaz_shape, blhaz_rate)
    ## LIKELIHOOD OF CURRENT DATA
    for ( i in 1:n[s] ) {
      linpred[s,i] <- x[s,i] * beta[s]
      for ( j in 1:Ji[s,i] )
        event[s,i,j] ~ dpois( blhaz[s, j] * exp( linpred[s,i] + logrisktime[s,i,j] ) )
    }
  }
})

## Non-informative prior code
noninf <- nimbleCode({
  ## LOOP THROUGH NUMBER OF STRATA
  for ( s in 1:S ) {
    ## PRIORS FOR EXCHANGEABLE PARAMETERS
    ## Prior for exchangeable treatment effects
    beta[s] ~ dnorm(beta_mean, sd = beta_sd)
    for ( j in 1:J[s] )
      blhaz[s,j] ~ dgamma(blhaz_shape, blhaz_rate)
    ## LIKELIHOOD OF CURRENT DATA
    for ( i in 1:n[s] ) {
      linpred[s,i] <- x[s,i] * beta[s]
      for ( j in 1:Ji[s,i] )
        event[s,i,j] ~ dpois( blhaz[s, j] * exp( linpred[s,i] + logrisktime[s,i,j] ) )
    }
  }
})