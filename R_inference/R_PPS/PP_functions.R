#---------------------------------------------------------------
# FUNCTIONS for the posterioir predictions of Labiodentals
#---------------------------------------------------------------
library(tidyr)


get_mean_LB<- function(yrep, phoible.area.sel){
  HGvsLB <- tibble(Area=unique(phoible.area.sel$Area))
  for (i in 1:nrow(yrep)){
    phoible.sim <- bind_cols(Area=phoible.area.sel$Area, Labiodentals.sim=yrep[i,],
                             has.Labiodentals.sim=(yrep[i,]!=0) )
    
    sim.lb <- phoible.sim %>%  dplyr::group_by(Area) %>% dplyr::summarize(Mean.lb=mean(has.Labiodentals.sim))
    sim.lb <- sim.lb[match(HGvsLB$Area, sim.lb$Area),]
    # rescale means as normalized
    rescaled <- (sim.lb$Mean.lb - mean(sim.lb$Mean.lb) )/sd(sim.lb$Mean.lb)
    sim.lb$Mean.lb <- rescaled
    ###
    HGvsLB <- bind_cols(HGvsLB, Mean.lb=sim.lb$Mean.lb)
    phoible.sim <- c()
  }
  
  return(HGvsLB)
}

####
do_regression <- function(mean.tb, Hg.mean, nsim){
  M1.coeff <- tibble(Intercept=rep(NA, nsim), slope.Hg=rep(NA, nsim))
  
  for (i in 1:nsim){
    dt.temp <- bind_cols(Lb=mean.tb[[1+i]], Hg=Hg.mean$mean)
    M1 <- glm(Lb ~ Hg, data = dt.temp, family = gaussian)
    M1.coeff[i,1] <-M1$coefficients[1]
    M1.coeff[i,2] <-M1$coefficients[2]
  }
  
  return(M1.coeff)
}


do_regression_BI <- function(mean.tb, Hg.mean, nsim){
  M1.coeff <- c()
  Hg.mean <-Hg.mean[match(mean.tb$Area, Hg.mean$Area),]
  
  for (i in 1:nsim){
    dt.temp <- bind_cols(Lb=mean.tb[[1+i]], Hg=Hg.mean$mean)
    
    M1 <- glm(Lb ~ Hg, data = dt.temp, family = gaussian)
    
    M1 <-stan_glm(Lb ~ Hg, data = dt.temp, family = gaussian, chains = 1, cores = 1)
    slope <- as.matrix(M1)[,2]
    
    M1.coeff <-c(M1.coeff, slope)
    
  }
  
  M1.coeff <-tibble(slope=M1.coeff)
  return(M1.coeff)
}

# Empirical distribution table for K-L divergence
get_dist_table <- function(breaks = seq(-10, 10, by=0.5), sample ){
  duration.cut = cut(sample, breaks, right=FALSE)
  duration.freq = table(duration.cut)
  duration.dist <- duration.freq/sum(duration.freq)
  return(duration.dist)
}
##
#------ FUNCTIONs for Error and Plotting

reshape_MeanTable_list <- function(across.area, MeanTable){
  reshaped.MeanTable <- c()
  for (i in seq_along(MeanTable)){
    out <- reshape_MeanTable(across.area, MeanTable[[i]] )
    out <- out %>% mutate(Dataset=names(MeanTable[i]))
    reshaped.MeanTable <- bind_rows(reshaped.MeanTable, out)
  }
  return(reshaped.MeanTable)
}


reshape_MeanTable <- function(across.area, mean.tb){
  mean.tb.resh <-across.area[,c(1,9)]
  mean.tb.resh <- left_join(mean.tb.resh, mean.tb)
  mean.tb.resh <-mean.tb.resh  %>% gather('sim', 'val', -Area, -mean.hg.emp.norm)
  return(mean.tb.resh)
}

#----- Mean sq error
mean_sq_erro <- function(across.area, mean.tb, type=c('common')){
  xx <- left_join(across.area[,c(1,8)], mean.tb)
  error <- apply(xx[,-c(1,2)], 2,  function(y) (y-xx$mean.lb.emp.norm)^2 ) %>%
    apply(., 1, function(yy) mean(yy))
  
  if (type=='common')
    return(mean(error))
  
  if (type=='individual')
    return( mutate(across.area[,1], error=error) )
}
#---------------------------------------------------------------------------