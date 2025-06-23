library(dplyr)
library(ggplot2)

load("data/mle.Rdata")

grid = expand.grid(p = c(0.0, 0.5, 1.0),
                   n = c(25, 50, 100),
                   q = c(-0.4, -0.2, 0, 0.2, 0.4))

# ## Save summary dataset
# res.all = data.frame()
# for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
#   res.all =
#     lapply(1:nrow(grid), function(hist) {
#       f = paste0('results/', method, '/res_summary_', hist, '.rds')
#       if (file.exists(f))
#         readRDS(f) %>%
#         mutate(q=grid$q[hist], p=grid$p[hist], n=grid$n[hist], method=method)
#     }) %>%
#     bind_rows() %>%
#     rbind(res.all, .)
# }
# res.all = res.all %>%
#   group_by(variable, q, p, n) %>%
#   mutate(rel_abs_bias = abs_bias / abs_bias[match('bhm', method)],
#          rel_mse = mse / mse[match('bhm', method)],
#          rel_ci_cov = ci_cov / 0.95,
#          rel_ci_width = ci_width / ci_width[match('bhm', method)])
# saveRDS(res.all,'results/sim_summary.rds')

res.all = readRDS('results/sim_summary.rds')

##########################################################
res.all = res.all %>%
  # filter(!(method=='bhm' & n!=100)) %>%
  mutate(method = case_when(method=='pp' ~ 'NPP + BHM',
                            method=='sleap' ~ 'SLEAP + BHM',
                            method=='pp2' ~ 'NPP',
                            method=='sleap2' ~ 'SLEAP',
                            TRUE ~ method))

# PLOT 1
# res.sub = res.all %>% filter(variable=='beta_mean')
# plot_name = 'plot_all_n'
# n.linetypes = rep(1, 4)
# n.linewidths = c(0.7, 1.4, 2.1, 2.8)
# n.alphas = c(0.4, 0.5, 0.6, 0.7)

# PLOT 2
res.sub = res.all %>% filter(variable=='beta_mean')
plot_name = 'plot'
n.linetypes = 3:1
n.linewidths = rep(0.8, 3)
n.alphas = rep(0.7, 3)

# PLOT 3
# res.sub = res.all %>% filter(variable=='beta[1]',
#                              n %in% c(25, 50, 100))
# plot_name = 'plot_beta_1'
# n.linetypes = 3:1
# n.linewidths = rep(0.8, 4)
# n.alphas = rep(0.7, 4)

# PLOT 4
# res.sub = res.all %>% filter(variable=='beta[7]',
#                              n %in% c(25, 50, 100))
# plot_name = 'plot_beta_7'
# n.linetypes = 3:1
# n.linewidths = rep(0.8, 4)
# n.alphas = rep(0.7, 4)

# PLOT 5
# res.sub = res.all %>% filter(variable=='beta[4]',
#                              n %in% c(25, 50, 100))
# plot_name = 'plot_beta_4'
# n.linetypes = 3:1
# n.linewidths = rep(0.8, 4)
# n.alphas = rep(0.7, 4)
##########################################################

## Make main plots
plots = list()
for (i in 1:2) { # for each p
  p.this = c(0, 0.5)[i]
  dat = res.sub %>% filter(p==p.this,
                           method!='sleap2')
  
  for (j in 1:3) { # for each result type
    type = c('rel_mse', 'rel_ci_cov', 'rel_ci_width')[j]
    type_lab = c('Rel. MSE', 'Rel. CI Coverage', 'Rel. CI Width')[j]
    
    plot = 
      ggplot(dat,
           aes(q, 
               eval(parse(text=type)),
               color=toupper(method),
               linetype=factor(n))) +
      geom_line(aes(linewidth=factor(n), alpha=factor(n))) +
      theme_minimal() +
      scale_linetype_manual('n per strata', values=n.linetypes) +
      scale_linewidth_manual('n per strata', values=n.linewidths) +
      scale_alpha_manual('n per strata', values=n.alphas) +
      scale_color_manual('method', values=c('NPP'='#66C2A5', 'BHM'='#FC8D62', 'NPP + BHM'='#1985A1', 'SLEAP'="#DDCC77", 'RMAPP.50'='#A6A670', 'RMAPP.75'='#86836D', 'SLEAP + BHM'="#882255", 'NONINF'="#BB4430")) +
      ggtitle(paste0('Proportion Exchangeable = ', p.this)) +
      labs(y = type_lab) +
      theme(plot.title = element_text(hjust = 0.5, size=10, face='bold'),
            legend.key.width = unit(2,"cm"),
            plot.margin = margin(0.3,0.3,0.3,0.3, "cm"),
            legend.position='none')
    
    plots[[(i-1)*3+j]] = ggplotGrob(plot)
  }
}

## Extract legend only - n
plot = 
  ggplot(dat,
         aes(q, 
             eval(parse(text=type)),
             linetype=factor(n))) +
  geom_line(aes(linewidth=factor(n), alpha=factor(n))) +
  theme_minimal() +
  scale_linetype_manual('n per strata', values=n.linetypes) +
  scale_linewidth_manual('n per strata', values=n.linewidths) +
  scale_alpha_manual('n per strata', values=n.alphas) +
  ggtitle(paste0('Proportion Exchangeable = ', p.this)) +
  labs(y = type_lab) +
  theme(plot.title = element_text(hjust = 0.5, size=10, face='bold'),
        legend.key.width = unit(2,"cm"),
        plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))
plots[[7]] = cowplot::get_legend(plot)

## Extract legend only - method
plot = 
  ggplot(dat,
         aes(q, 
             eval(parse(text=type)),
             color=toupper(method))) +
  geom_line(linewidth=0.8, alpha=0.7) +
  theme_minimal() +
  scale_color_manual('method', values=c('NPP'='#66C2A5', 'BHM'='#FC8D62', 'NPP + BHM'='#1985A1', 'SLEAP'="#DDCC77", 'RMAPP.50'='#A6A670', 'RMAPP.75'='#86836D', 'SLEAP + BHM'="#882255", 'NONINF'="#BB4430")) +
  ggtitle(paste0('Proportion Exchangeable = ', p.this)) +
  labs(y = type_lab) +
  theme(plot.title = element_text(hjust = 0.5, size=10, face='bold'),
        legend.key.width = unit(2,"cm"),
        plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))
plots[[8]] = cowplot::get_legend(plot)


## Combine plots and save
jpeg(paste0('results/', plot_name, '.jpeg'), width = 7500, height = 5000, res=600)
grid.arrange(grobs=plots, layout_matrix=rbind(c(1,1,2,2,3,3),
                                              c(1,1,2,2,3,3),
                                              c(4,4,5,5,6,6),
                                              c(4,4,5,5,6,6),
                                              c(7,8,NULL,NULL,NULL,NULL)))
dev.off()


## Save Bias plot separately
plots2 = list()
for (i in 1:2) { # for each p
  p.this = c(0, 0.5)[i]
  dat = res.sub %>% filter(p==p.this)
  
  type = 'bias'
  type_lab = 'Bias'
  
  plot = 
    ggplot(dat,
           aes(q, 
               eval(parse(text=type)),
               color=toupper(method),
               linetype=factor(n))) +
    geom_line(aes(linewidth=factor(n), alpha=factor(n))) +
    theme_minimal() +
    scale_linetype_manual('n per strata', values=n.linetypes) +
    scale_linewidth_manual('n per strata', values=n.linewidths) +
    scale_alpha_manual('n per strata', values=n.alphas) +
    scale_color_manual('method', values=c('NPP'='#66C2A5', 'BHM'='#FC8D62', 'NPP + BHM'='#1985A1', 'SLEAP'="#DDCC77", 'RMAPP.50'='#A6A670', 'RMAPP.75'='#86836D', 'SLEAP + BHM'="#882255", 'NONINF'="#BB4430")) +
    ggtitle(paste0('Proportion Exchangeable = ', p.this)) +
    labs(y = type_lab) +
    theme(plot.title = element_text(hjust = 0.5, size=10, face='bold'),
          legend.key.width = unit(2,"cm"),
          plot.margin = margin(0.3,0.3,0.3,0.3, "cm"),
          legend.position='none')
  
  plots2[[i]] = ggplotGrob(plot)
}
plots2[[3]] = plots[[7]] # legend - n
plots2[[4]] = plots[[8]] # legend - method

jpeg(paste0('results/', plot_name, '_bias.jpeg'), width = 5000, height = 3000, res=600)
grid.arrange(grobs=plots2, layout_matrix=rbind(c(1,1,2,2),
                                               c(1,1,2,2),
                                               c(3,4,NULL,NULL)))
dev.off()
