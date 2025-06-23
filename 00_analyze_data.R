library(dplyr)
library(survival)

## Get data
c = bind_rows(read.csv('data/surv_sanofi.csv') %>% rename_all(tolower) %>%
                rename(subjid=rsubjid) %>% mutate(histdata=1, kras=NA),
              read.csv('data/surv_amgen.csv') %>% rename_all(tolower) %>%
                mutate(histdata=0) %>% filter(!kras %in% c('Failed', '', 'Mutant'))) %>% # wild-type (trt efficacious) group only
  select(!c(starts_with(c("pfs",  "pdeath")))) %>%
  mutate(os=os/12)

## Create strata key-table
strata.table = c %>% group_by(ecog2, bevacizumab, oxaliplatin) %>%
  summarize(n0=sum(histdata==1), n=sum(histdata==0)) %>% arrange(desc(n0)) %>%
  ungroup() %>% mutate(s=row_number(), j=ifelse(ceiling(n/20)<=5, ceiling(n/20), 5))

c = c %>% left_join(strata.table%>%select(!c(n0,n,j)), by=c("ecog2", "bevacizumab", "oxaliplatin"))

S = nrow(strata.table)
J = strata.table$j
n = strata.table$n
n0 = strata.table$n0
trt.prop = mean(c$trt[c$histdata==0])

## MLE within strata in pooled historical + current data
trteff = rep(NA, S)
blshape = rep(NA, S)
scale = rep(NA, S)
event.rate.trt = rep(NA, S)
event.rate.ctrl = rep(NA, S)
for (s in 1:S) {
  d = c[c$s==s,]
  
  event.rate.trt[s] = mean(d$death[d$trt==1])
  event.rate.ctrl[s] = mean(d$death[d$trt==0])
  
  # Fit Weibull
  fit = survreg(Surv(os, death)~trt, data=d, dist='weibull')
  blshape[s] = 1 / fit$scale
  scale[s] = exp(coef(fit)[1]) ^(-blshape[s])
  trteff[s] = -coef(fit)[2]
}

trteff = c(-0.3341, -0.2602, -0.1622, -0.1420, -0.1023, 0.0299, 0.2702)

save(list=c('S', 'J', 'n', 'n0', 'strata.table', 
            'trt.prop', 'blshape', 'scale', 'trteff',
            'event.rate.trt', 'event.rate.ctrl'),
     file="data/mle.Rdata")
