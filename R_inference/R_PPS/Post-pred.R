###########################
#
# Posterior prediction for across group variation in the probability of labiodentals and hunter-gatherers
# Comment 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################
library("rstanarm")
library("bayesplot")

#---- set directory
# set working directory to 'Tarasov_Uyeda_Supplementary' then
setwd("./R_inference/R_PPS")

# Read necessary objects and data ------------

# handy functions
source('PP_functions.R')

# data obtained at CM_plot_pars_estimates.R, include parameter estimated for H and L
Par.tb.hg <- readRDS(file = '../R_CM/output/Par.tb.hg.rds') # Par.tb.hg
Par.tb.lb <- readRDS(file = '../R_CM/output/Par.tb.lb.rds') # Par.tb.lb

# phoible.area (phoible_distance_dat.csv) is the same as the original phoible but includes distance from Africa
# this distance is obtained as described in PP_get_geo-dist.R
phoible.area<-read.csv("output/phoible_distance_dat.csv", header=T, stringsAsFactors = F, sep=',')
phoible.area<-as_tibble(phoible.area)
phoible.area

#------ make datatset 
library("plyr")
phoible.area <-phoible.area %>% dplyr::mutate(Afr.HG= replace(rep(0,n()), which(Area=='Africa' & Subsistence.TOM=='HG'), 1 )  )
phoible.area <- phoible.area %>% dplyr::mutate(NC.HG= replace(rep(0,n()), which(Area=='N-C Asia' & Subsistence.TOM=='HG'), 1 )  )

# calculate the residual phonemes NonLb.NonFr.phonemes as (Nonlabiodental.phonemes-Nonlabiodental.Fricatives)
phoible.area <- phoible.area %>% dplyr::mutate(NonLb.NonFr.phonemes=Nonlabiodental.phonemes-Nonlabiodental.Fricatives)
phoible.area <- phoible.area %>% dplyr::mutate(Subsistence.TOM.num= mapvalues(Subsistence.TOM, c('AGR', 'HG'), c(0,1) ) %>% as.numeric() )



###########################################################################
# TOM DATATSET ------------------------------------------------------------
###########################################################################

# filter data for PP simulations using TOM --------------
phoible.area.sel <- phoible.area %>% dplyr::select(Labiodentals,  Nonlabiodental.Fricatives, Area, Dist, 
                                                   Subsistence.TOM, Subsistence.TOM.num, NonLb.NonFr.phonemes, has.labiodentals)

# remove 'Africa' and 'N-C Asia' since they support dependent model in Causal Models inference
phoible.area.sel <- phoible.area.sel %>% filter(Area!='Africa' & Area!='N-C Asia')
phoible.area.sel <-phoible.area.sel%>% na.omit()

#--------------------------------------------------------

# Create across.area object that summurizes empirical probabilities of H and L
#  and those inffered using causal models
ml <- Par.tb.lb %>% filter(Dataset=='Tom')
ml <-ml[1:10,]
mh <- Par.tb.hg %>% filter(Dataset=='Tom' & Model=='Ind')
mh <-mh[1:10,]

across.area <- tibble(Area=ml$Area, mean.lb=ml$mean, mean.hg=mh$mean)
# noprmalize data
across.area <- across.area %>% mutate(mean.lb.norm= (mean.lb-mean(mean.lb))/sd(mean.lb) )
across.area <- across.area %>% mutate(mean.hg.norm= (mean.hg-mean(mean.hg))/sd(mean.hg) )

# remove Africa and N-C Asia
across.area <- across.area[-c(1,10),]

#-- Empirical means 
# Empirical means of the parameters correspond to their max Ln estimates
emp.lb <- phoible.area.sel %>%  dplyr::group_by(Area) %>% dplyr::summarize(mean.lb.emp=mean(has.labiodentals))
emp.hg <-  phoible.area.sel %>%  dplyr::group_by(Area) %>% dplyr::summarize(mean.hg.emp=mean(Subsistence.TOM.num))
across.area <- left_join(across.area, emp.lb)
across.area <- left_join(across.area, emp.hg)
# normalize
across.area <- across.area %>% mutate(mean.lb.emp.norm= (mean.lb.emp-mean(mean.lb.emp))/sd(mean.lb.emp) )
across.area <- across.area %>% mutate(mean.hg.emp.norm= (mean.hg.emp-mean(mean.hg.emp))/sd(mean.hg.emp) )

# check how BI estimates differ from the empirical ones
across.area$mean.lb-across.area$mean.lb.emp
across.area$mean.hg-across.area$mean.hg.emp


# Regression on observed data ---------------------------------------------
library("ggridges")

# Regression on the parameters (H and G) estimated using binomial causal models
biR0.est <- stan_glm(mean.lb.norm ~ 1, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000)
biR1.est <- stan_glm(mean.lb.norm ~ 0+mean.hg.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000)
biR2.est <- stan_glm(mean.lb.norm ~ mean.hg.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000)
summary(biR0.est)
summary(biR1.est)
summary(biR2.est)

# Regression on the empirical probailities of H and G

biR0.emp.est <- stan_glm(mean.lb.emp.norm ~ 1, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r0.csv"))
biR1.emp.est <- stan_glm(mean.lb.emp.norm ~ 0+mean.hg.emp.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r1.csv"))
biR2.emp.est <- stan_glm(mean.lb.emp.norm ~ mean.hg.emp.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r2.csv"))

# Marginal ML
ml0 <- bridge_sampler(biR0.emp.est)
#Bridge sampling estimate of the log marginal likelihood: -15.61533
ml1 <- bridge_sampler(biR1.emp.est)
#Bridge sampling estimate of the log marginal likelihood: -12.3594
ml2 <- bridge_sampler(biR2.emp.est)
#Bridge sampling estimate of the log marginal likelihood:  -15.84426
bayes_factor(ml1, ml0, log=T)
# BF 3.25592

summary(biR0.emp.est)
summary(biR1.emp.est) 
# Estimates:
#                    mean   sd    2.5%   25%   50%   75%   97.5%
# mean.hg.emp.norm  -0.6    0.3  -1.3   -0.8  -0.6  -0.4   0.0  
# sigma              0.9    0.3   0.5    0.7   0.8   1.0   1.5  
# mean_PPD           0.0    0.3  -0.6   -0.2   0.0   0.2   0.6  
# log-posterior    -11.8    1.1 -14.9  -12.2 -11.5 -11.0 -10.7 

summary(biR2.emp.est)

# -----------------------------------
# PLOTS -----------------------------

# Compare slopes for the empirical and estimated probabilities of H and L

ESvsEMP <- bind_rows(
  tibble(vals=as.matrix(biR1.est)[,1], model='Estimated'),
  tibble(vals=as.matrix(biR1.emp.est)[,1], model='Empirical')
)
ggplot(ESvsEMP , aes(x=vals, y=model)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7)

# Compare slopes for empirical biR1.emp.est and biR2.emp.est models
biR1.emp.slope <- as_tibble(biR1.emp.est)
biR <- bind_rows(
  tibble(vals=as.matrix(biR1.emp.est)[,1], model='R1'),
  tibble(vals=as.matrix(biR2.emp.est)[,2], model='R2')
)
ggplot(biR, aes(x=vals, y=model)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7)


###### Plot regression using empirical H and L

Pp <- ggplot(across.area) + 
  #ggplot(across.area) + 
  aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm) + 
  geom_abline(aes(intercept = 0, slope = mean.hg.emp.norm), 
              data = biR1.emp.slope, color = "grey60", 
              alpha = .03) +
  geom_abline(intercept = 0, 
              slope = mean(biR1.emp.slope$mean.hg.emp.norm), 
              size = 1, color = "#3366FF") +
  geom_point()  +
  geom_text(aes(label=Area),hjust=.4, vjust=2, size=3) +
  labs(x = 'Pr of hunter-gatherers', y = 'Pr of labiodental phonemes') +
  scale_x_continuous(limits = c(-1.3, 2.3))+
  scale_y_continuous(limits = c(-1.3, 2))

#ggsave(paste0('output/Pr_HGvsLB_tom', '.pdf'), Pp, width = 100, height = 110, units = c('mm') )
ggsave(paste0('output/Pr_HGvsLB_tom', '.png'), Pp, width = 100, height = 110, units = c('mm') )


#-----------------------------------------
# TOM, Posterior Predictions: SIMULATIONS
#-----------------------------------------

nsim=2000
# hg mean to use in the simulated regressions
Hg.mean <- across.area %>% select(Area, mean=mean.hg.emp.norm)

# objects to store the info
MeanTable.Tom <- vector(length = 6, mode='list')
names(MeanTable.Tom) <- c('Distance', 'Nonlabiodental Fricatives', 'Nonlabiodental Phonemes', 'Subsistence', 'Area', 'Intercept')

# --- Dist
LbDist <- stan_glm(Labiodentals ~ Dist, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[1]] <- mean.tb


# --- Nonlabiodental.Fricatives
LbDist <- stan_glm(Labiodentals ~ Nonlabiodental.Fricatives, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[2]] <- mean.tb

# --- NonLb.NonFr.phonemes
LbDist <- stan_glm(Labiodentals ~ NonLb.NonFr.phonemes, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[3]] <- mean.tb

# --- Subsistence.TOM
LbDist <- stan_glm(Labiodentals ~ Subsistence.TOM, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
summary(LbDist)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[4]] <- mean.tb

# --- Area
LbDist <- stan_glm(Labiodentals ~ Area, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[5]] <- mean.tb

# --- Randomly simulated Labiodentals without any predictors (only Intercept)
LbDist <- stan_glm(Labiodentals ~ 1, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.Tom[[6]] <- mean.tb


#  Calculate Mean Squared Error -------
Error.Tom <- lapply(MeanTable.Tom, function(z) mean_sq_erro(across.area, z, type='common') )
Error.Tom <-c(Dataset='TOM', Error.Tom)
Error.Tom <-as_tibble(Error.Tom)
Error.Tom 
xtable(Error.Tom, type = "latex", digits = 3, caption ='Mean Sq Error; TOM')


# --- PLOT
reshaped.MeanTable.Tom <- reshape_MeanTable_list(across.area,  MeanTable.Tom)


#theme_get()
#theme_set(theme_grey())

Pp <- ggplot(across.area) + 
  geom_point(data=reshaped.MeanTable.Tom, aes(x=mean.hg.emp.norm, y=val, color=Area), alpha=.2, size=.7, stroke=0) +
  geom_point(aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm), alpha=1, size=4)+
  geom_point(aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm, color=Area), alpha=1, size=3)+
  
  geom_text(aes(label=Area, y = mean.lb.emp.norm, x = mean.hg.emp.norm),hjust=.4, vjust=2, size=3) +
  labs(x = 'Pr of hunter-gatherers', y = 'Pr of labiodental phonemes') +
  scale_x_continuous(limits = c(-1.3, 2.3))+
  scale_y_continuous(limits = c(-1.5, 2)) +
  facet_wrap(~ Dataset, nrow = 3, ncol = 2) +
  theme(legend.position = "none")
Pp

ggsave(paste0('output/PP_Lb-per-area_tom', '.png'), Pp, width = 300, height = 200, units = c('mm') )
##################################################################################################







###########################################################################
# AUTOTYP DATATSET -------------------------------------------------------
###########################################################################

## phoible data for the  simulations using AUTOTYP

phoible.area.sel <- phoible.area %>% dplyr::select(Labiodentals,  Nonlabiodental.Fricatives, Area, Dist, 
                                                   Subsistence.AUTOTYP, NonLb.NonFr.phonemes, has.labiodentals)
phoible.area.sel <-  phoible.area.sel %>% dplyr::mutate(Subsistence.AUTOTYP.num= mapvalues(Subsistence.AUTOTYP, c('AGR', 'HG'), c(0,1) ) %>% as.numeric() )
phoible.area.sel <- phoible.area.sel %>% filter(Area!='Africa' & Area!='N-C Asia')
phoible.area.sel <- phoible.area.sel %>%na.omit()

# -- make data
ml <- Par.tb.lb %>% filter(Dataset=='Auto')
ml <-ml[1:10,]
mh <- Par.tb.hg %>% filter(Dataset=='Auto' & Model=='Ind')
mh <-mh[1:10,]

across.area <- tibble(Area=ml$Area, mean.lb=ml$mean, mean.hg=mh$mean)
# noprmalize data
across.area <- across.area %>% mutate(mean.lb.norm= (mean.lb-mean(mean.lb))/sd(mean.lb) )
across.area <- across.area %>% mutate(mean.hg.norm= (mean.hg-mean(mean.hg))/sd(mean.hg) )

# remove Africa and N-C Asia
across.area <- across.area[-c(1,10),]

#-- Empirical means 
# Empirical means of the parameters correspond to their max Ln estimates
emp.lb <- phoible.area.sel %>%  dplyr::group_by(Area) %>% dplyr::summarize(mean.lb.emp=mean(has.labiodentals))
emp.hg <-  phoible.area.sel %>%  dplyr::group_by(Area) %>% dplyr::summarize(mean.hg.emp=mean(Subsistence.AUTOTYP.num))
across.area <- left_join(across.area, emp.lb)
across.area <- left_join(across.area, emp.hg)
# normalize
across.area <- across.area %>% mutate(mean.lb.emp.norm= (mean.lb.emp-mean(mean.lb.emp))/sd(mean.lb.emp) )
across.area <- across.area %>% mutate(mean.hg.emp.norm= (mean.hg.emp-mean(mean.hg.emp))/sd(mean.hg.emp) )

# check how BI estimates differ from the empirical ones
across.area$mean.lb-across.area$mean.lb.emp
across.area$mean.hg-across.area$mean.hg.emp


#------------------ Regression on observed data
library("ggridges")

# Regression on the empirical estimates of H and LB
biR0.emp.est <- stan_glm(mean.lb.emp.norm ~ 1, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r0.csv"))
biR1.emp.est <- stan_glm(mean.lb.emp.norm ~ 0+mean.hg.emp.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r1.csv"))
biR2.emp.est <- stan_glm(mean.lb.emp.norm ~ mean.hg.emp.norm, data = across.area, family = gaussian, chains = 2, cores = 1, iter = 5000,
                         diagnostic_file = file.path(tempdir(), "df_r2.csv"))
# Marginal ML
ml0 <- bridge_sampler(biR0.emp.est)
#Bridge sampling estimate of the log marginal likelihood: -15.62577
ml1 <- bridge_sampler(biR1.emp.est)
#Bridge sampling estimate of the log marginal likelihood: -13.15717
ml2 <- bridge_sampler(biR2.emp.est)
#Bridge sampling estimate of the log marginal likelihood: -16.54143
bayes_factor(ml1, ml0, log=T)
# BF 2.46860

summary(biR0.emp.est)
summary(biR1.emp.est)
# Estimates:
#   mean   sd    2.5%   25%   50%   75%   97.5%
# mean.hg.emp.norm  -0.5    0.4  -1.3   -0.7  -0.5  -0.3   0.3  
# sigma              1.0    0.3   0.6    0.8   0.9   1.1   1.7  
# mean_PPD           0.0    0.4  -0.7   -0.2   0.0   0.2   0.7  
# log-posterior    -12.7    1.2 -15.7  -13.1 -12.4 -11.9 -11.6

summary(biR2.emp.est)


biR1.emp.slope <- as_tibble(biR1.emp.est)
####

####### PLOTS

# compare slopes in biR1.emp.est and biR2.emp.est models
biR <- bind_rows(
  tibble(vals=as.matrix(biR1.emp.est)[,1], model='R1'),
  tibble(vals=as.matrix(biR2.emp.est)[,2], model='R2')
)
ggplot(biR, aes(x=vals, y=model)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7)


# Plot regression
Pp <- ggplot(across.area) + 
  aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = 0, slope = mean.hg.emp.norm), 
              data = biR1.emp.slope, color = "grey60", 
              alpha = .03) +
  # Plot the median values in blue
  geom_abline(intercept = 0, 
              slope = mean(biR1.emp.slope$mean.hg.emp.norm), 
              size = 1, color = "#3366FF") +
  geom_point()  +
  geom_text(aes(label=Area),hjust=.4, vjust=2, size=3) +
  labs(x = 'mean hunter-gatherers', y = 'mean labiodental phonemes') +
  scale_x_continuous(limits = c(-1.3, 2.3))+
  scale_y_continuous(limits = c(-1.3, 2))

#ggsave(paste0('output/Pr_HGvsLB__autotyp', '.pdf'), Pp, width = 210, height = 100, units = c('mm') )
ggsave(paste0('output/Pr_HGvsLB_autotyp', '.png'), Pp, width = 210, height = 100, units = c('mm') )


#-----------------------------------------
# AUTOTYP, Posterior Predictions: SIMULATIONS
#-----------------------------------------

nsim=2000
# hg mean to use in the simulated regressions
Hg.mean <- across.area %>% select(Area, mean=mean.hg.emp.norm)

# objects to store the info
MeanTable.AUTO <- vector(length = 6, mode='list')
names(MeanTable.AUTO) <- c('Distance', 'Nonlabiodental Fricatives', 'Nonlabiodental Phonemes', 'Subsistence', 'Area', 'Intercept')

# --- Dist
LbDist <- stan_glm(Labiodentals ~ Dist, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[1]] <- mean.tb


# --- Nonlabiodental.Fricatives
LbDist <- stan_glm(Labiodentals ~ Nonlabiodental.Fricatives, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[2]] <- mean.tb

# --- NonLb.NonFr.phonemes
LbDist <- stan_glm(Labiodentals ~ NonLb.NonFr.phonemes, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[3]] <- mean.tb

# --- Subsistence.TOM
LbDist <- stan_glm(Labiodentals ~ Subsistence.TOM, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
summary(LbDist)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[4]] <- mean.tb

# --- Area
LbDist <- stan_glm(Labiodentals ~ Area, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[5]] <- mean.tb

# --- Randomly simulated Labiodentals without any predictors
LbDist <- stan_glm(Labiodentals ~ 1, data = phoible.area.sel, family = poisson, chains = 2, cores = 1, iter = 5000)
yrep<- posterior_predict(LbDist, draws = nsim)
mean.tb <- get_mean_LB(yrep, phoible.area.sel)
MeanTable.AUTO[[6]] <- mean.tb



# ------- Sq Error
Error.AUTO <- lapply(MeanTable.AUTO, function(z) mean_sq_erro(across.area, z, type='common') )
Error.AUTO <-c(Dataset='AUTO', Error.AUTO)
Error.AUTO <-as_tibble(Error.AUTO)

# combine TOM and AUTO
Error <- bind_rows(Error.Tom, Error.AUTO)
Error
# Dataset Distance `Nonlabiodental Fricatives` `Nonlabiodental Phonemes` Subsistence  Area Intercept
# <chr>      <dbl>                       <dbl>                     <dbl>       <dbl> <dbl>     <dbl>
#   1 TOM        0.437                        1.04                      1.44       0.861 0.143      1.72
# 2 AUTO       0.824                        1.51                      1.66       1.61  0.523      1.75

xtable(Error, type = "latex", digits = 3, caption ='Mean Squared Error')


# --- PLOT
reshaped.MeanTable.AUTO <- reshape_MeanTable_list(across.area,  MeanTable.AUTO)

theme_get()
theme_set(theme_grey())

Pp <- ggplot(across.area) + 
  geom_point(data=reshaped.MeanTable.AUTO, aes(x=mean.hg.emp.norm, y=val, color=Area), alpha=.2, size=.7, stroke=0) +
  geom_point(aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm), alpha=1, size=4)+
  geom_point(aes(y = mean.lb.emp.norm, x = mean.hg.emp.norm, color=Area), alpha=1, size=3)+
  
  geom_text(aes(label=Area, y = mean.lb.emp.norm, x = mean.hg.emp.norm),hjust=.4, vjust=2, size=3) +
  labs(x = 'Pr of hunter-gatherers', y = 'Pr of labiodental phonemes') +
  scale_x_continuous(limits = c(-1.7, 2))+
  scale_y_continuous(limits = c(-1.3, 2)) +
  facet_wrap(~ Dataset, nrow = 3, ncol = 2) +
  theme(legend.position = "none")
Pp

ggsave(paste0('output/PP_Lb-per-area_auto', '.png'), Pp, width = 300, height = 200, units = c('mm') )



