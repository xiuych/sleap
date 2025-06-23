
#' Weibull PH density
#' @param x where to evaluate density
#' @param a treatment indicator. Must be an vector of 0s and 1s of length 1 or length(x)
#' @param shape shape parameter for Weibull
#' @param scale scale parameter for Weibull
#' @param beta log hazard ratio
#' @param log whether to return log of density function
#' @return evaluation of Weibull PH density at each x
dweibph <- function(x, a, shape, scale, beta, log = FALSE) {
  if ( !(length(a) %in% c(1, length(x) ) ) )
    stop('a must be of length 1 or the same dimension as x')
  if ( !( all(a %in% c(0, 1) ) ) )
    stop('a must be equal to 0 or 1')
  logx        <- log(x)
  logscale    <- log(scale) + beta * a
  logf        <- log(shape) + logscale + (shape - 1) * logx - exp(logscale + shape * logx)
  logf[x < 0] <- -Inf
  if ( log )
    return(logf)
  return(exp(logf))
}

#' Weibull PH survival function
#' @param q quantile
#' @param a treatment indicator. Must be an vector of 0s and 1s of length 1 or length(x)
#' @param shape shape parameter for Weibull
#' @param scale scale parameter for Weibull
#' @param beta log hazard ratio
#' @param lower.tail whether to return distribution (TRUE) or survival (FALSE) function
#' @param log.p whether to return log of distribution function
#' @return evaluation of Weibull PH distribution function at each x
pweibph <- function(q, a, shape, scale, beta, lower.tail = FALSE, log.p = FALSE) {
  if ( !(length(a) %in% c(1, length(q) ) ) )
    stop('a must be of length 1 or the same dimension as q')
  if ( !( all(a %in% c(0, 1) ) ) )
    stop('a must be equal to 0 or 1')
  logq         <- log(q)
  logscale     <- log(scale) + beta * a
  logS         <- -exp(logscale + shape * logq)
  logS[q <= 0] <- 0
  if ( lower.tail )
    logS <- logS + log(expm1(-logS))   ## = log(1 - S(t))
  if ( log.p )
    return(logS)
  return(exp(logS))
}

#' Weibull PH random number generation (via inverse CDF)
#' @param n number of variables to generate
#' @param a treatment indicator. Must be an vector of 0s and 1s of length 1 or n
#' @param shape shape parameter for Weibull distribution
#' @param scale scale parameter for Weibull distribution
#' @param beta log hazard ratio
#' @return vector of random variables
rweibph <- function(n, a, shape, scale, beta) {
  if ( !(length(a) %in% c(1, n ) ) )
    stop('a must be of length 1 or n')
  if ( !( all(a %in% c(0, 1) ) ) )
    stop('a must be equal to 0 or 1')
  logscale <- log(scale) + beta * a
  Z        <- rexp(n)
  exp( (log(Z) - logscale) / shape )
}

#' Find optimal censoring rate for desired proportion of event times
#' @param lambda Rate parameter for Exponential censoring time
#' @param shape Weibull shape parameter
#' @param beta log hazard ratio
#' @param a treatment indicator. Must be a scalar equal to 0 or 1
#' @param prop.event desired proportion of observed event times. Should be between 0, 1.
#' @return optimal value of lambda
lambdaopt <- function(a, shape, scale, beta, prop.event) {
  optfun <- function(lambda) {
    ## Compute probability C > T
    int <- integrate(function(x) {
      exp( 
        pexp(x, lambda, lower.tail = F, log.p = T) 
        + dweibph(x, a, shape, scale, beta, log = TRUE) 
      )
    }
    , 0, Inf )
    Pevent <- int$value
    abs(Pevent - prop.event)
  }
  optimize(optfun, c(0, 100))
}


#' Generate data set
#' @param n number of subjects
#' @param prob.trt probability of receiving treatment (treatment assignment is not random)
#' @param shape shape parameter for Weibull distribution
#' @param scale scale parameter for Weibull distribution
#' @param beta regression coefficient (log hazard ratio)
#' @param cens.rate.trt rate of censoring among treated individuals. If NULL, will compute censoring rate based on prob.event.trt
#' @param cens.rate.ctrl rate of censoring among controls. If NULL, will compute censoring rate based on prob.event.ctrl
#' @param prob.event.trt if `prob.rate.ctrl` is NULL, will compute censoring rate to match `prob.event.trt`
#' @param prob.event.ctrl if `prob.rate.ctrl` is NULL, will compute censoring rate to match `prob.event.trt`
#' @return complete data set
gendata <- function(
    n, prob.trt, shape, scale, beta, cens.rate.trt = NULL, cens.rate.ctrl = NULL
    , prop.event.trt = NULL, prop.event.ctrl = NULL
) {
  ## Generate treatment
  ntrt <- ceiling(n * prob.trt)
  a    <- c(rep(0, n - ntrt), rep(1, ntrt))
  # a <- rbinom(n, 1, prob.trt)
  ## Generate Weibull event time
  t <- rweibph(n, a, shape, scale, beta)
  ## Generate censoring time
  if ( is.null(cens.rate.trt) )
    cens.rate.trt <- lambdaopt(a = 1, shape, scale, beta, prop.event.trt)$minimum
  if ( is.null(cens.rate.ctrl) )
    cens.rate.ctrl <- lambdaopt(a = 0, shape, scale, beta, prop.event.ctrl)$minimum
  c <- numeric(n)
  c[a==0] <- rexp(sum(a==0), cens.rate.ctrl)
  c[a==1] <- rexp(sum(a==1), cens.rate.trt)
  ## Create observed time
  y     <- pmin(t, c)
  event <- t <= c
  ## Create data set
  dat <- data.frame(y = y, event = event, a = a, t = t, c = c)
  attr(dat, 'params') <- c(
    'prob.trt' = prob.trt, 'shape' = shape, 'scale' = scale
    , 'beta' = beta, 'cens.rate.trt' = cens.rate.trt, 'cens.rate.ctrl' = cens.rate.ctrl
  )
  return(dat)
}

#' Negative log likelihood of Weibull PH model
#' @param parms current value of parmaeters. Takes for (log(shape), log(scale), beta)
#' @param y observed event times
#' @param event event indicators
#' @param ctrlonly logical indicating whether to fit model with only controls (i.e., no beta)
#' @return negative of log likelihood
negloglik <- function(parms, y, event, a, ctrlonly = FALSE) {
  shape  <- exp( parms[1] )
  scale  <- exp( parms[2] )
  if (ctrlonly) {
    beta <- 0
  } else {
    beta <- parms[3]
  }
  indx   <- event==1
  loglik <- numeric(length(y))
  loglik[indx]  <- dweibph(y[indx], a[indx], shape, scale, beta, log = TRUE)
  loglik[!indx] <- pweibph(y[!indx], a[!indx], shape, scale, beta, lower.tail = FALSE, log.p = TRUE)
  -sum(loglik)
}
#' Gradient of the negative log likelihood for a Weibull PH model
gnegloglik <- function(parms, y, event, a, ctrlonly = FALSE) {
  shape  <- exp( parms[1] )
  scale  <- exp( parms[2] )
  res    <- numeric(length(parms))
  if (ctrlonly) {
    beta <- 0
  } else {
    beta <- parms[3]
  }
  indx    <- event==1
  logy    <- log(y)
  nevents    <- sum(event)
  ntrtevents <- sum(a * event)
  Q1  <- exp(beta * a + shape * logy)
  Q2  <- Q1 * logy
  Q3  <- a * Q1      
  
  # res[1] <- nevents / shape + sum(event * logy) - scale * sum(exp(beta * a) * y^shape * log(y))
  res[1] <- nevents / shape + sum(event * logy) - scale * sum(Q2)
  res[2] <- nevents / scale - sum(Q1)
  if(!ctrlonly)
    res[3] <- sum(a * event) - scale * sum(Q3)
  ## Chain rule
  res[1] <- res[1] * shape
  res[2] <- res[2] * scale
  -res
}
