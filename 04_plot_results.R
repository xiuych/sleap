library(shiny)
library(dplyr)
library(ggplot2)
library(survival)
library(ggfortify)

load("data/mle.Rdata")

grid = expand.grid(p = c(0.0, 0.5, 1.0),
                   n = c(25, 50, 100),
                   q = c(-0.4, -0.2, 0, 0.2, 0.4))

res.all = data.frame()
for (method in c('sleap', 'pp', 'bhm', 'noninf')) {
  res.all = 
    lapply(1:nrow(grid), function(hist) {
      f = paste0('results/', method, '/res_summary_', hist, '.rds')
      if (file.exists(f))
        readRDS(f) %>%
        mutate(q=grid$q[hist], p=grid$p[hist], n=grid$n[hist], method=method)
    }) %>%
    bind_rows() %>%
    rbind(res.all, .)
}

ui <- fluidPage(
  tabsetPanel(tabPanel('Method Comparison',
                       column(12, align = "right", actionButton('update', 'Update (slurm jobs will be submitted via you)',
                                                                style='margin-top:10px;')),
                       p(paste0('beta_truth = (', paste(as.character(round(trteff, 2)), collapse=', '), ')')),
                       selectInput(inputId = "var",
                                   label = "Variable",
                                   choices = unique(res.all$variable),
                                   selected = 'beta_mean'),
                       selectInput(inputId = "p",
                                   label = "p",
                                   choices = unique(res.all$p)),
                       selectInput(inputId = "n",
                                   label = "n",
                                   choices = unique(res.all$n)),
                       selectInput(inputId = "measure",
                                   label = "Measure",
                                   choices = c('abs_bias_pct', 'mse', 'ci_cov', 'ci_width'),
                                   selected = 'ci_cov'),
                       plotOutput('plot_methods')),
              tabPanel('LEAP Exchangeability',
                       selectInput(inputId = 'gamma',
                                   label = 'Variable',
                                   choices=grep('^gamma', unique(res.all$variable), value=T)),
                       selectInput(inputId = "n2",
                                   label = "n",
                                   choices = unique(res.all$n)),
                       plotOutput('plot_exch')),
              tabPanel('Historical Data Info',
                       selectInput(inputId = 'hist',
                                   label = 'Scenario',
                                   choices = 1:nrow(grid)),
                       textOutput('histinfo'),
                       plotOutput('survplot')))
)

server <- function(input, output) {
  output$plot_methods <- renderPlot({
    ggplot(res.all %>% filter(variable==input$var, p==input$p, n==input$n),
           aes(q, 
               eval(parse(text=input$measure)),
               color=toupper(method))) +
      geom_line(linewidth=0.8, alpha=.7) +
      ggtitle(paste0('median nsims=', median(filter(res.all, variable==input$var, p==input$p)$nsims))) +
      theme_minimal() +
      scale_color_manual('method', values=c('PP2'='#66C2A5', 'BHM'='#FC8D62', 'PP'='#1985A1', 'SLEAP2'="#DDCC77", 'RMAPP.50'='#A6A670', 'RMAPP.75'='#86836D', 'SLEAP'="#882255", 'NONINF'="#BB4430")) +
      labs(y = c('Absolute Bias (%)', 'MSE', 'CI Coverage')[which(c('abs_bias_pct', 'mse', 'ci_cov')==input$measure)])
  })
  
  observeEvent(input$update, {
    system(paste0('cd ', getwd()), intern=T)
    terminalout = system(paste0('sbatch ', getwd(), '/03_compile_sims.sh'), intern=T)
    showNotification(terminalout)
    })
  
  output$plot_exch <- renderPlot({
    ggplot(res.all %>% filter(variable==input$gamma, n==input$n2),
           aes(q,
               eval(parse(text=c('mean_postmean'))),
               color=toupper(p))) +
      geom_line(linewidth=0.8, alpha=.6) +
      theme_minimal() +
      scale_color_manual('p', values=c('#BB4430', '#7EBDC2', '#F3DFA2')) +
      labs(y = c('mean_postmean +- sd_postmean')) +
      geom_errorbar(aes(ymin=eval(parse(text=c('mean_postmean')))-eval(parse(text=c('sd_postmean'))),
                        ymax=eval(parse(text=c('mean_postmean')))+eval(parse(text=c('sd_postmean')))), width=.1)
  })
  
  output$histinfo <- renderText({
    a = readRDS(paste0('data/historical/', 'histdata_', input$hist, '.rds'))
    paste0('Baseline hazard errors: ', round(a$err, 3),
           ', p = ', toString(a$p),
           ', q = ', a$q)
  })
  
  output$survplot <- renderPlot({
    a = readRDS(paste0('data/historical/', 'histdata_', input$hist, '.rds'))$data
    autoplot(survfit(Surv(y, death) ~ s, data = a))
  })
}

shinyApp(ui,server)
