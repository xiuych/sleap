Packages <- c("ggplot2", "survival", "ks", 
              #"survminer", 
              "MASS", "hesim", 
              #"rjags", 
              #"psborrow", 
              "cmdstanr",
              "pracma")
lapply(Packages, library, character.only = TRUE)

idx_samp = 1
idx_scale = 1
idx_n=1
idx_a0_n   <- 1
cut.time.mod <- c(0, sim.data %>% filter(nu == 1) %>% pull(t) %>% median) 
n.seg          <- length(cut.time.mod) + 1
interval.names <- paste0("interval", 1:length(cut.time.mod))
cov.names      <- NULL

scale_list <- NA
T1E_vec <- NA
X3_vec <- NA
cen_vec <- NA
x_t_vec <- NA
h_n_vec <- NA

## Below code adapted from Kwiatkowski et al 01B-config-output-containers.R

################################################################################
### Part 4: Output args ########################################################
################################################################################

a <- c("HR", "MSE", "lower", "upper", "CI_width", "CI", "power"#, "se", "BIC", cov.names
)
b <- c("a0", "0", "0.25", "0.5", "0.75", "1", "cox")
c <- c("com_cauchy", "com_a0")
d <- c("a0_v1")

# pqList <- c("p1_q0", "p3_q0", "p5_q0", "p7_q0")
# pqList <- paste0("p", format(seq(1, 5, by = 0.1), nsmall = 1))
pqList <- paste0("p", format(3, nsmall = 1))
# tList <- paste("t", format(c(0), nsmall = 3), sep = "_")
# tList <- paste("t", format(c(0.950, 0.955, 0.960, 0.965, 0.970, 0.975,
#                              0.990, 0.991, 0.993, 0.995, 0.997, 0.999, 0), nsmall = 3), sep = "_")
tList <- paste("t", format(c(0.990, 0), nsmall = 3), sep = "_")
pqList <- paste(expand.grid(pqList, tList)$Var1, expand.grid(pqList, tList)$Var2, sep = "_")
nList  <- c(b, pqList)
ppList <- c(nList[1:6], pqList)

names <- c(paste(expand.grid(a,b)$Var1, expand.grid(a,b)$Var2, sep = "_"),
           paste(expand.grid(a,c)$Var1, expand.grid(a,c)$Var2, sep = "_"),
           paste(expand.grid(a,d)$Var1, expand.grid(a,d)$Var2, sep = "_"),
           "SE_a0", "a0_logHR", "bias_a0", "mean_a0", 
           c(paste0("a0_wt_int", 1:(n.seg - 1))),
           c(paste(expand.grid(c("a0_mean", paste0("a0_wt_int", 1:(n.seg - 1))), c("q0.1", "q0.25", "q0.5", "q0.75", "q0.9"))$Var1,
                   expand.grid(c("a0_mean", paste0("a0_wt_int", 1:(n.seg - 1))), c("q0.1", "q0.25", "q0.5", "q0.75", "q0.9"))$Var2, sep = "_")),
           "mean_a0_r", "a0_mean_r_q0.1", "a0_mean_r_q0.25", "a0_mean_r_q0.5", "a0_mean_r_q0.75", "a0_mean_r_q0.9",
           "nu_trt", "nu_ctrl", "nu_h", 
           "a0_wt_obs", "a0_wt_cen",
           paste(expand.grid(a,pqList)$Var1, expand.grid(a,pqList)$Var2, sep = "_"),
           "scale", "T1E", "X3", "cen", "x_t", "h_n")

################################################################################
### Part 5) Set containers #####################################################
################################################################################

outer_table <- matrix(NA, nrow = length(scale_list), ncol = length(names))
colnames(outer_table) <- names
outermost <- matrix(NA, nrow = idx_samp, ncol = length(names))
colnames(outermost) <- names
# outer_a0 <- matrix(NA, nrow = (n.seg - 1) * max(h_n_vec), ncol = idx_samp)

################################################################################
### Part 6: Transformation functions ###########################################
################################################################################

# rotate3 <- function(x, p){
#    ((2 * (x - 0.5)) ^ p + 1) / 2
# }

rotate3 <- function(x, p){
  2 ^ (p - 1) * sign(x - 0.5) * abs(x - 0.5) ^ p + 0.5
}

g <- function(x, q){
  1 / (1 + exp(-(200*(x - q))))
}




## Below code based on Kwiatkowski et al, 02A1-fit-initial-models.R

################################################################################
################################################################################
# EVENTS
################################################################################
################################################################################

# 2023-11-06 equal events per interval
cut.time.mod <- c(0, sim.data %>% filter(nu == 1) %>% pull(t) %>% median) 
inner.t        <- cut.time.mod[-1] # asssuming n.seg >= 3
diff           <- cut.time.mod[-1] - cut.time.mod[-length(cut.time.mod)] # asssuming n.seg >= 3
# cut.time.mod
sim.split          <- survSplit(sim.data, 
                                cut = cut.time.mod,
                                # cut = cut.time, 
                                end = "t", 
                                start = "t0", 
                                event = "nu", 
                                episode = "interval")
sim.split$interval <- factor(sim.split$interval - 1)
sim.split$expo     <- sim.split$t - sim.split$t0

### Loop 1) (idx_outer) fit models for historical events/cen and RCT events ####
## FORMULAS FOR EVENTS ALWAYS INVOLVE COVARIATES
formula2 <- formula(paste("nu ~ treat +", paste(cov.names, collapse=" + "), " interval + offset(log(expo)) - 1"))
formula1 <- formula(paste("nu ~ treat +", paste(cov.names, collapse=" + "), " offset(log(expo))"))
formula2_no_trt <- formula(paste("nu ~ ", paste(cov.names, collapse=" + "), " interval + offset(log(expo)) - 1"))
formula1_no_trt <- formula(paste("nu ~ ", paste(cov.names, collapse=" + "), " offset(log(expo))"))
if (n.seg >= 3){
  # Fit model for historical events
  sim.poisson <- glm(formula2_no_trt, 
                     family = poisson(link = "log"), 
                     data = sim.split[sim.split$external == 1, ])
  # Fit model for RCT events
  sim.poisson.lambda <- glm(formula2,
                            family = poisson(link = "log"), 
                            data = sim.split[sim.split$external == 0, ])
} else {
  # Fit model for historical events
  sim.poisson <- glm(formula1_no_trt, 
                     family = poisson(link = "log"), 
                     data = sim.split[sim.split$external == 1, ])
  names(sim.poisson$coefficients)[names(sim.poisson$coefficients) == "(Intercept)"] <- "interval1"
  # Fit model for RCT events
  sim.poisson.lambda <- glm(formula1,
                            family = poisson(link = "log"), 
                            data = sim.split[sim.split$external == 0, ])
  names(sim.poisson.lambda$coefficients)[names(sim.poisson.lambda$coefficients) == "(Intercept)"] <- "interval1"
}

################################################################################
################################################################################
# INSTANCES OF CENSORING
################################################################################
################################################################################

# 2023-11-19 equal events per interval
cut.time.cen <- c(0, sim.data %>% filter(nu == 0 & external == 1) %>% pull(t) %>% median)
inner.t.cen  <- cut.time.cen[-1] # asssuming n.seg >= 3
diff.cen     <- cut.time.cen[-1] - cut.time.cen[-length(cut.time.cen)] # asssuming n.seg >= 3
# cut.time.cen
sim.split.cen          <- survSplit(sim.data, 
                                    cut = cut.time.cen,
                                    # cut = cut.time, 
                                    end = "t", 
                                    start = "t0", 
                                    event = "nu", 
                                    episode = "interval")
sim.split.cen$interval <- factor(sim.split.cen$interval - 1)
sim.split.cen$expo     <- sim.split.cen$t - sim.split.cen$t0

sim.split$interval.cen <- 1 + as.numeric(inner.t.cen > inner.t) # 2023-11-19 only works with n.seg == 3

# head(sim.split)
# head(sim.split.cen)

### Loop 1) (idx_outer) fit models for historical events/cen and RCT events ####
if (n.seg >= 3){
  # Fit model for historical censoring
  sim.poisson.theta   <- glm(1 - nu ~ interval + offset(log(expo)) - 1,
                             family = poisson(link = "log"),
                             data = sim.split.cen[sim.split.cen$external == 1, ])
} else {
  # Fit model for historical censoring
  sim.poisson.theta   <- glm(1 - nu ~ offset(log(expo)),
                             family = poisson(link = "log"),
                             data = sim.split.cen[sim.split.cen$external == 1, ])
  names(sim.poisson.theta$coefficients)[names(sim.poisson.theta$coefficients) == "(Intercept)"] <- "interval1"
}



## Below code based on Kwiatkowski et al, 02-full-sims.R


for(idx_outer in 1){
  set.seed(23 * as.integer(mod(idx_scale, 800) * idx_samp + idx_outer))  # 2020-05-19
  
  ### Loop 2) (idx) # reps per design (i.e. num of imputations per exposure time interval in external controls)
  innermost  <- matrix(NA, nrow = idx_n, ncol = length(names))
  colnames(innermost) <- names
  a0.outer <- matrix(NA, nrow = nrow(sim.split[sim.split$external == 1, ]), ncol = idx_n)
  for(idx in 1:idx_n){
    lambda <- mvrnorm(1, sim.poisson$coefficients[c(cov.names, interval.names)], 
                      vcov(sim.poisson)[c(cov.names, interval.names), 
                                        c(cov.names, interval.names)])
    theta <- mvrnorm(1, sim.poisson.theta$coefficients[interval.names],
                     vcov(sim.poisson.theta)[interval.names, interval.names])
    
    ### Loop 2) (idx) Add imputation for multiple intervals in historical controls
    sim.split.idx <- sim.split # 2021-08-03
    if (n.seg >= 3){ 
      for(i in which(sim.split.idx$external == 1)){
        for(j in 1:length(inner.t)){
          if (sim.split.idx$nu[i] == 0 & sim.split.idx$interval[i] == j & sim.split.idx$expo[i] == diff[j]){
            temp.t <- rexp(1, rate = exp(as.matrix(sim.split.idx[i, c(cov.names)]) %*%
                                           as.matrix(lambda[c(cov.names)]) + 
                                           lambda[paste0("interval", j)]))
            # temp.c  <- rexp(1, rate = exp(theta[paste0("interval", j)]))
            temp.c  <- rexp(1, rate = exp(theta[paste0("interval", sim.split.idx$interval.cen[i])]))
            temp.y  <- pmin(temp.t, temp.c)
            temp.nu <- as.numeric(temp.y == temp.t)
            sim.split.idx$expo[i] <- diff[j] + temp.y
            sim.split.idx$nu[i]   <- temp.nu # NOT EVEN USED AGAIN 2023-11-19
            
            if (is.na(sim.split.idx$expo[i])){
              print(i)
              print(j)
            }
            
          }
        }
      }
    }
    
    # subset to historical only
    sim.split.h <- sim.split.idx[sim.split.idx$external == 1, ]
    box.p.log <- matrix(NA, ncol = idx_a0_n, nrow = nrow(sim.split.h))
    
    ### Loop 3) (idx_a0) Repetitions per a0 (i.e. number of compatibility scores computed per exposure time interval)
    for (idx_a0 in 1:idx_a0_n){
      lambda <- mvrnorm(1, sim.poisson.lambda$coefficients[c(cov.names, interval.names)], 
                        vcov(sim.poisson.lambda)[c(cov.names, interval.names), 
                                                 c(cov.names, interval.names)])
      theta <- mvrnorm(1, sim.poisson.theta$coefficients[interval.names],
                       vcov(sim.poisson.theta)[interval.names, 
                                               interval.names])
      # debugging code - prevent theta from being positive
      for (i in 1:length(theta)){
        if (theta[i] > 0) theta[i] <- min(theta)
      }
      
      ### Loop 3) (idx_a0) generate y.h.pred incorporating variability in estimating lambda and assess compatibility, use interval + expo
      for(i in 1:nrow(sim.split.h)){
        j <- c(sim.split.h[i, "interval"])
        
        if (exp(as.matrix(sim.split.h[i, c(cov.names)]) %*%
                as.matrix(lambda[c(cov.names)]) +
                lambda[paste0("interval", j)]) == 0) {## Coco edit
          t.h.pred=Inf
        } else {
          t.h.pred <- rexp(1E4, rate = exp(as.matrix(sim.split.h[i, c(cov.names)]) %*%
                                             as.matrix(lambda[c(cov.names)]) +
                                             lambda[paste0("interval", j)]))
        }
        # c.h.pred <- rexp(1E4, rate = exp(theta[paste0("interval", j)]))
        if (j == 1){ # 2023-11-19 only works if n.seg == 3
          if (inner.t.cen >= inner.t){
            c.h.pred <- rexp(1E4, rate = exp(theta[paste0("interval",  j)]))
          } else {
            c.h.pred <- rpwexp(1E4, time = cut.time.cen, rate = exp(theta))
          }
        }
        if (j == 2){ # 2023-11-19 only works if n.seg == 3
          if (inner.t >= inner.t.cen){
            c.h.pred <- rexp(1E4, rate = exp(theta[paste0("interval",  j)]))
          } else {
            c.h.pred <- rpwexp(1E4, time = cut.time.cen - cut.time.mod, rate = exp(theta))
          }
        }
        
        y.h.pred <- pmin(t.h.pred, c.h.pred)
        # fhat_outer         <- kde(sqrt(y.h.pred))
        # fhat_predict_outer <- predict(fhat_outer, x = sqrt(y.h.pred))
        # box.p.log[i, idx_a0] <- mean(fhat_predict_outer <= predict(fhat_outer, x = sqrt(sim.split.h$expo[i])))
        fhat_outer         <- kde(log(y.h.pred))
        fhat_predict_outer <- predict(fhat_outer, x = log(y.h.pred))
        box.p.log[i, idx_a0] <- mean(fhat_predict_outer <= predict(fhat_outer, x = log(sim.split.h$expo[i])))
      }
    } # ends for (idx_a0 in 1:idx_a0_n)
    
    
    ### Loop 2) (idx) Compute case-specific weights
    a0 <- rowMeans(box.p.log)
    
    # debugging code - prevent element of a0 from being identically zero
    a0[a0 == 0] <- 1E-6
    
    a0.outer[, idx] <- a0
    sim.split.h$a0 <- a0
    # Compute average weights a0
    innermost[idx, "mean_a0"] <- mean(a0)
    innermost[idx, c("a0_mean_q0.1", "a0_mean_q0.25", "a0_mean_q0.5", "a0_mean_q0.75", "a0_mean_q0.9")] <- quantile(a0, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    innermost[idx, "mean_a0_r"] <- mean(rotate3(a0, p = 1))
    innermost[idx, c("a0_mean_r_q0.1", "a0_mean_r_q0.25", "a0_mean_r_q0.5", "a0_mean_r_q0.75", "a0_mean_r_q0.9")] <- quantile(rotate3(a0, p = 1), probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    innermost[idx, "a0_wt_obs"] <- mean(sim.split.h[sim.split.h$nu == 1, "a0"])
    innermost[idx, "a0_wt_cen"] <- mean(sim.split.h[sim.split.h$nu == 0, "a0"])
    # Compute average weights a0 by interval
    for (j1 in 1:(n.seg - 1)){
      innermost[idx, paste0("a0_wt_int", j1)] <- mean(sim.split.h[sim.split.h[,"interval"] == j1, "a0"])
      innermost[idx, paste0("a0_wt_int", j1, c("_q0.1", "_q0.25", "_q0.5", "_q0.75", "_q0.9"))] <- quantile(sim.split.h[sim.split.h[,"interval"] == j1, "a0"], probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    }
    
  } # ends for(idx in 1:idx_n) # repetitions per design (i.e. number of imputations per exposure time interval in external controls)
  
  
  ### Loop 1) (idx_outer) For each idx_scale and idx_sample, summarize data from idx_n repeats
  innermost <- data.frame(innermost)
  outermost[idx_outer, ] <- colMeans(innermost)
  
  ### Loop 1) (idx_outer) Final analysis using fixed power prior weights & adaptive weights
  for (i in 1:length(ppList)){
    
    # create wts vector
    if (ppList[i] == "a0"){
      wts <- c(rep(1, times = sum(sim.split$external == 0)), rowMeans(a0.outer))
    } else if (substring(ppList[i], 1, 1) == "p"){
      wts <- c(rep(1, times = sum(sim.split$external == 0)), 
               rotate3(x = rowMeans(a0.outer), 
                       p = as.numeric(substring(ppList[i], 2, 4))) *
                 g(x = 2 * mean(a0.outer), 
                   q = as.numeric(substring(ppList[i], 8, 12))))
    } else {
      wts <- c(rep(1, times = sum(sim.split$external == 0)), rep(as.numeric(ppList[i]), sum(sim.split$external == 1)))
    }
    
    if (n.seg >= 3){
      sim.poisson <- glm(formula2,
                         family = poisson(link = "log"),
                         data = sim.split,
                         weights = wts)
    } else {
      sim.poisson <- glm(formula1,
                         family = poisson(link = "log"),
                         data = sim.split,
                         weights = wts)
    }
    
    outermost[idx_outer, paste0(c("HR_", "lower_", "upper_"), ppList[i])] <-
      c(exp(coef(sim.poisson))["treat"],
        exp(coef(sim.poisson)["treat"])*exp(-qnorm(0.975)*summary(sim.poisson)$coefficients["treat", 2]),
        exp(coef(sim.poisson)["treat"])*exp(qnorm(0.975)*summary(sim.poisson)$coefficients["treat", 2]))
    if (ppList[i] == "a0"){
      outermost[idx_outer, "a0_logHR"] <- coef(sim.poisson)["treat"]
      outermost[idx_outer, "SE_a0"] <- summary(sim.poisson)$coefficients["treat", 2]
      temp = matrix(c(exp(coef(sim.poisson)["treat"]) * summary(sim.poisson)$coefficients["treat", 2],
                      confint(sim.poisson)["treat",]), nrow=1)
      colnames(temp) <- c("HR_SE_a0", "lower_logHR_a0", "upper_logHR_a0")
      outermost = outermost %>% cbind(temp)
    }
  }
  
  
  ### Loop 1) (idx_outer) Final analysis using cox model
  cox_formula <- formula(paste("Surv(t, nu) ~ treat ", paste(cov.names, collapse=" + ")))
  cox <- coxph(cox_formula, data = sim.data[sim.data$external == 0, ])
  # outermost[idx_outer, c("HR_cox", "lower_cox", "upper_cox", paste0(cov.names, "_cox"))] <-
  #   c(exp(coef(cox))["treat"],
  #     exp(coef(cox)["treat"])*exp(-qnorm(0.975)*summary(cox)$coefficients["treat", "se(coef)"]),
  #     exp(coef(cox)["treat"])*exp(qnorm(0.975)*summary(cox)$coefficients["treat", "se(coef)"]),
  #     summary(cox)$coefficients[cov.names, "coef"])
  outermost[idx_outer, c("HR_cox", "lower_cox", "upper_cox")] <-
    c(exp(coef(cox))["treat"],
      exp(coef(cox)["treat"])*exp(-qnorm(0.975)*summary(cox)$coefficients["treat", "se(coef)"]),
      exp(coef(cox)["treat"])*exp(qnorm(0.975)*summary(cox)$coefficients["treat", "se(coef)"]))
  
  
  
  
  ### Loop 1) (idx_outer) Compute MSE, HR, credible interval width, and power for each analysis method
  for (i in 1:length(nList)){
    # outermost[, paste0("MSE_", nList[i])]      <- 
    #   (outermost[, paste0("HR_", nList[i])] - exp(mu2["treat"])) ^ 2
    outermost[, paste0("CI_width_", nList[i])] <- 
      (outermost[, paste0("upper_", nList[i])] - outermost[, paste0("lower_", nList[i])])
    # outermost[, paste0("CI_", nList[i])]       <- 
    #   (outermost[, paste0("lower_", nList[i])] < exp(mu2["treat"]) & 
    #      exp(mu2["treat"]) < outermost[, paste0("upper_", nList[i])])
    outermost[, paste0("power_", nList[i])]    <- 
      ((outermost[, paste0("lower_", nList[i])] < 1 & outermost[, paste0("upper_", nList[i])] < 1) 
       # | (outermost[, paste0("lower_", nList[i])] > 1 & outermost[, paste0("upper_", nList[i])] > 1)
      )
  }
  
  outermost[idx_outer, "nu_trt"]  <- sum(sim.data$nu[sim.data$treat == 1 & sim.data$external == 0])
  outermost[idx_outer, "nu_ctrl"] <- sum(sim.data$nu[sim.data$treat == 0 & sim.data$external == 0])
  outermost[idx_outer, "nu_h"]    <- sum(sim.data$nu[sim.data$treat == 0 & sim.data$external == 1])
  
  # scale_list[idx_scale]
  # T1E_vec[idx_scale]
  # X3_vec[idx_scale]
  # cen_vec[idx_scale]
  # x_t_vec[idx_scale]
  # h_n_vec[idx_scale]
  
  outermost[idx_outer, "scale"] <- scale_list[idx_scale]
  outermost[idx_outer, "T1E"]   <- T1E_vec[idx_scale] 
  outermost[idx_outer, "X3"]    <- X3_vec[idx_scale] 
  outermost[idx_outer, "cen"]   <- cen_vec[idx_scale] 
  outermost[idx_outer, "x_t"]   <- x_t_vec[idx_scale]
  outermost[idx_outer, "h_n"]   <- h_n_vec[idx_scale] 
  
  # outer_a0[1:length(a0), idx_outer] <- a0
  # print(idx_outer)

} # ends for(idx_outer in 1:idx_samp)


res_s = outermost[,c('a0_logHR', 'SE_a0', 'lower_logHR_a0', 'upper_logHR_a0',
                     'HR_a0', 'HR_SE_a0', 'lower_a0', 'upper_a0')]
