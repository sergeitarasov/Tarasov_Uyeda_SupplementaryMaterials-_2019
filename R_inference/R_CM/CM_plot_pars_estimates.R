#==================================================
#
# Plot paramter estimates for Causal Models
#
# Tarasov, Uyeda 2019
#=================================================

# needs the followimg objects contiainng BI estimates
#fit.tom.dep
#fit.tom.ind
#fit.auto.dep
#fit.auto.ind

#==================
#
# Plot Parameters and save them to files
#
#==================
library(bayesplot)
library("gridExtra")
setwd("./output")


# TOM  Dependent models

#== TOM specify model info
obj <- fit.tom.dep
model.tag <- 'Dependent-model_Tom-dataset'
#==

param.plots <- vector(13, mode="list")
for (i in 1:length(obj)){
  
  plot.mcmc <- as.matrix(obj[[i]])
  
  param.plots[[i]] <- mcmc_intervals(
    plot.mcmc, 
    pars = c("theta_S",  "theta_L1", "theta_L2"),
    prob = 0.95, # 90% intervals
    prob_outer = 1, # 100%
    point_est = "mean"
  ) +
    ggplot2::labs(
      title = names(obj[i])
      #,subtitle = model.tag
    ) +
    ggplot2::theme(plot.title = element_text(size=10), axis.text.y=element_text(size=8) ) +
    scale_y_discrete( labels=c('HG', 'LB(HG)', 'LB(AGR)') )
}

many.plots <- marrangeGrob(param.plots, nrow=5, ncol=3,  top = NULL )
#ggsave(paste0('Pars_', model.tag, '.pdf'), many.plots, width = 210, height = 297, units = c('mm') )
ggsave(paste0('Pars_', model.tag, '.png'), many.plots, width = 210, height = 297, units = c('mm') )


#==================
# AUTO Dependent models

#== AUTO specify model info
obj <- fit.auto.dep
model.tag <- 'Dependent-model_AUTO-dataset'
#==

param.plots <- vector(13, mode="list")
for (i in 1:length(obj)){
  
  plot.mcmc <- as.matrix(obj[[i]])
  
  param.plots[[i]] <- mcmc_intervals(
    plot.mcmc, 
    pars = c("theta_S",  "theta_L1", "theta_L2"),
    prob = 0.95, # 90% intervals
    prob_outer = 1, # 100%
    point_est = "mean"
  ) +
    ggplot2::labs(
      title = names(obj[i])
      #,subtitle = model.tag
    ) +
    ggplot2::theme(plot.title = element_text(size=10), axis.text.y=element_text(size=8) ) +
    scale_y_discrete( labels=c('HG', 'LB(HG)', 'LB(AGR)') )
}

many.plots <- marrangeGrob(param.plots, nrow=5, ncol=3,  top = NULL )
#ggsave(paste0('Pars_', model.tag, '.pdf'), many.plots, width = 210, height = 297, units = c('mm') )
ggsave(paste0('Pars_', model.tag, '.png'), many.plots, width = 210, height = 297, units = c('mm') )





#==================
# TOM  INDependent models
#fit.tom.ind[[1]]@model_pars

#== TOM specify model info
obj <- fit.tom.ind
model.tag <- 'Independent-model_Tom-dataset'
#==

param.plots <- vector(13, mode="list")
for (i in 1:length(obj)){
  
  plot.mcmc <- as.matrix(obj[[i]])
  
  param.plots[[i]] <- mcmc_intervals(
    plot.mcmc, 
    pars = c("theta_S",  "theta_L1"),
    prob = 0.95, # 90% intervals
    prob_outer = 1, # 100%
    point_est = "mean"
  ) +
    ggplot2::labs(
      title = names(obj[i])
      #,subtitle = model.tag
    ) +
    ggplot2::theme(plot.title = element_text(size=10), axis.text.y=element_text(size=8) ) +
    scale_y_discrete( labels=c('HG', 'LB') )
}

many.plots <- marrangeGrob(param.plots, nrow=5, ncol=3,  top = NULL )
#ggsave(paste0('Pars_', model.tag, '.pdf'), many.plots, width = 210, height = 297, units = c('mm') )
ggsave(paste0('Pars_', model.tag, '.png'), many.plots, width = 210, height = 297, units = c('mm') )

#==================
# AUTO INDependent models

#== AUTO specify model info
obj <- fit.auto.ind
model.tag <- 'Independent-model_AUTO-dataset'
#==

param.plots <- vector(13, mode="list")
for (i in 1:length(obj)){
  
  plot.mcmc <- as.matrix(obj[[i]])
  
  param.plots[[i]] <- mcmc_intervals(
    plot.mcmc, 
    pars = c("theta_S",  "theta_L1"),
    prob = 0.95, # 90% intervals
    prob_outer = 1, # 100%
    point_est = "mean"
  ) +
    ggplot2::labs(
      title = names(obj[i])
      #,subtitle = model.tag
    ) +
    ggplot2::theme(plot.title = element_text(size=10), axis.text.y=element_text(size=8) ) +
    scale_y_discrete( labels=c('HG', 'LB') )
}

many.plots <- marrangeGrob(param.plots, nrow=5, ncol=3,  top = NULL )
#ggsave(paste0('Pars_', model.tag, '.pdf'), many.plots, width = 210, height = 297, units = c('mm') )
ggsave(paste0('Pars_', model.tag, '.png'), many.plots, width = 210, height = 297, units = c('mm') )


#==================================================
#
# Make Tables
#
# 
#=================================================
#fit.tom.dep
#fit.tom.ind
#fit.auto.dep
#fit.auto.ind


#==================
# HG parameters

#== specify parameter
param <- 'HG'
#==

tb <- c()
tb.names <- c()
for (i in 1:13){
  
  block <- rbind(
    summary(fit.tom.ind[[i]])$summary[1,c(1,3,4,8)],
    summary(fit.tom.dep[[i]])$summary[1,c(1,3,4,8)],
    summary(fit.auto.ind[[i]])$summary[1,c(1,3,4,8)],
    summary(fit.auto.dep[[i]])$summary[1,c(1,3,4,8)]
  )
  
  tb <-rbind(tb, block)
  
  tb.names <-rbind(tb.names,
  cbind(names(fit.tom.ind[i]), c('Tom', 'Tom', 'Auto', 'Auto'), c('Ind', 'Dep', 'Ind', 'Dep') )
  #cbind(names(fit.tom.ind[i]), c('T', 'T', 'A', 'A'), c('I', 'D', 'I', 'D') )
  )
  #c('Ti', 'Td', 'Ai', 'Ad')
}
colnames(tb.names) <- c('Area', 'Dataset', 'Model')
#

Par.tb <- bind_cols(as_tibble(tb.names), as_tibble(tb) )
Par.tb.hg <- Par.tb

Pp <- ggplot(Par.tb, aes(mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, x=paste(Dataset, Model), y=mean, color=Area ) ) +
  theme(legend.position="none", plot.title = element_text(size=8), axis.text.x=element_text(size=8, angle=90),
        strip.text.x = element_text(size = 8, angle=90)  )+
  scale_x_discrete( labels=c('Ti', 'Td', 'Ai', 'Ad') )+
  labs(x = NULL, y=param)+
  facet_grid(.~Area ) 
  #facet_wrap(~Area, nrow=5)

#ggsave(paste0(param, '_plot', '.pdf'), Pp, width = 210, height = 100, units = c('mm') )
ggsave(paste0(param, '_plot', '.png'), Pp, width = 210, height = 100, units = c('mm') )



#==================
# LB parameters

#== specify parameter
param <- 'LB'
#==

tb <- c()
tb.names <- c()
for (i in 1:13){
  
  block <- rbind(
    summary(fit.tom.ind[[i]])$summary[2,c(1,3,4,8)],
    #summary(fit.tom.dep[[i]])$summary[1,c(1,3,4,8)],
    summary(fit.auto.ind[[i]])$summary[2,c(1,3,4,8)]
    #summary(fit.auto.dep[[i]])$summary[1,c(1,3,4,8)]
  )
  
  tb <-rbind(tb, block)
  
  tb.names <-rbind(tb.names,
                   cbind(names(fit.tom.ind[i]), c('Tom', 'Auto') )
                   #cbind(names(fit.tom.ind[i]), c('T', 'T', 'A', 'A'), c('I', 'D', 'I', 'D') )
  )
  #c('Ti', 'Td', 'Ai', 'Ad')
}
colnames(tb.names) <- c('Area', 'Dataset')
#

Par.tb <- bind_cols(as_tibble(tb.names), as_tibble(tb) )
Par.tb.lb <- Par.tb

Pp <- ggplot(Par.tb, aes(mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, x=Dataset, y=mean, color=Area ) ) +
  theme(legend.position="none", plot.title = element_text(size=8), axis.text.x=element_text(size=8, angle=90),
        strip.text.x = element_text(size = 8, angle=90)  )+
  scale_x_discrete( labels=c('A', 'T') )+
  labs(x = NULL, y=param)+
  facet_grid(.~Area ) 


#ggsave(paste0(param, '_plot', '.pdf'), Pp, width = 210, height = 100, units = c('mm') )
ggsave(paste0(param, '_plot', '.png'), Pp, width = 210, height = 100, units = c('mm') )


#==================
# Lb(HG) and Lb(AGR) parameters

#== specify parameter
param <- 'Lb(HG) and Lb(AGR)'
#==

tb <- c()
tb.names <- c()
for (i in 1:13){
  
  block <- rbind(
    summary(fit.tom.dep[[i]])$summary[2,c(1,3,4,8)],
    summary(fit.auto.dep[[i]])$summary[2,c(1,3,4,8)],
    summary(fit.tom.dep[[i]])$summary[3,c(1,3,4,8)],
    summary(fit.auto.dep[[i]])$summary[3,c(1,3,4,8)]
  )
  
  tb <-rbind(tb, block)
  
  tb.names <-rbind(tb.names,
                   cbind(names(fit.tom.ind[i]), c('Tom', 'Auto', 'Tom', 'Auto'), c('Lh', 'Lh', 'La', 'La') )
                   #cbind(names(fit.tom.ind[i]), c('T', 'T', 'A', 'A'), c('I', 'D', 'I', 'D') )
  )
  #c('Ti', 'Td', 'Ai', 'Ad')
}
colnames(tb.names) <- c('Area', 'Dataset', 'Lb_param')
#

Par.tb <- bind_cols(as_tibble(tb.names), as_tibble(tb) )

Pp <- ggplot(Par.tb, aes(mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`, x=paste(Dataset, Lb_param), y=mean, color=Area ) ) +
  theme(legend.position="none", plot.title = element_text(size=8), axis.text.x=element_text(size=8, angle=90),
        strip.text.x = element_text(size = 7, angle=90)  )+
  #scale_x_discrete( labels=c('A', 'T') )+
  labs(x = NULL, y=param)+
  facet_grid(.~Area ) 


#ggsave(paste0('Lbs_dep_plot', '.pdf'), Pp, width = 210, height = 100, units = c('mm') )
ggsave(paste0('Lbs_dep_plot', '.png'), Pp, width = 210, height = 100, units = c('mm') )


#==================
# The following two object are needed for poserior predictions
Par.tb.hg
Par.tb.lb
#saveRDS(Par.tb.hg, file = 'Par.tb.hg.rds')
#saveRDS(Par.tb.lb, file = 'Par.tb.lb.rds')
