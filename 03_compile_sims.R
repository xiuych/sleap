print('r script started')
library(dplyr)
library(purrr)

load("data/mle.Rdata")


grid = expand.grid(p = c(0.0, 0.5, 1.0),
                   n = c(25, 50, 100),
                   q = c(-0.4, -0.2, 0, 0.2, 0.4))

id=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

for (hist in 1:nrow(grid)) {
  if (!is.na(id) && hist!=id) next
  if (grid$p[hist] == 1) next # no need to run p=1 scenario
  
  p = grid$p[hist]
  q = grid$q[hist]
  
  for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
    
    
    res_all = paste0('results/raw/', method) %>% dir(paste0('^res_', hist, '_'), full.names=T) %>%
      map_dfr(readRDS) %>%
      mutate(truth = case_when(
        variable=='beta_mean' ~ -0.10,
        variable=='beta_sd' ~ 0.20,
        variable=='beta[1]' ~ trteff[1],
        variable=='beta[2]' ~ trteff[2],
        variable=='beta[3]' ~ trteff[3],
        variable=='beta[4]' ~ trteff[4],
        variable=='beta[5]' ~ trteff[5],
        variable=='beta[6]' ~ trteff[6],
        variable=='beta[7]' ~ trteff[7],
        grepl('gamma\\[\\d, 1]', variable) ~ p
      )) %>%
      mutate(coverage = `2.5%`<=truth & truth<=`97.5%`)
    
    print('ok')
    
    res_summary = res_all %>%
      group_by(variable) %>%
      summarize(
        'nsims'         = sum(!is.na(mean))
        , 'truth'         = mean(truth)
        , 'mean_postmean' = mean(mean, na.rm = T)
        , 'sd_postmean'   = sd(mean, na.rm = T)
        , 'bias'          = mean(mean-truth, na.rm=T)
        , 'bias_pct'      = mean((mean - truth)/truth*100, na.rm = T)
        , 'abs_bias'      = mean(abs(mean - truth), na.rm=T)
        , 'abs_bias_pct'  = mean(abs((mean - truth)/truth*100), na.rm=T)
        , 'mse'           = mean((mean - truth)^2, na.rm = T)
        , 'ci_cov'        = mean(coverage, na.rm = T )
        , 'ci_width'      = mean(`97.5%` - `2.5%`, na.rm = T )
        , .groups='drop'
      )
    
    saveRDS(res_summary, file=paste0('results/', method, '/res_summary_', hist, '.rds'))
  }
}
