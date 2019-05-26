###########################
# Functions for making stan data for Phoible dataset.
#
# Fitting dependent and independent causal models for the language dataset from
# 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################


make_stan_data <- function(dt){
  ## make data
  N=sum(dt$n) # total number of observation
  
  # number of HG languegaes
  N.hg <- ( dt %>% filter(Subsistence.TOM=='HG') %>% summarise(sum(n)) %>% c())
  N.hg <- N.hg$`sum(n)`
  
  # number of HG languages with labiodentals
  L.hg <- dt %>% filter(Subsistence.TOM=='HG', has.labiodentals==T) %>% c()
  L.hg <-L.hg$n
  
  # number of AGR languages with labiodentals
  L.agr <- dt %>% filter(Subsistence.TOM=='AGR', has.labiodentals==T) %>% c()
  L.agr <-L.agr$n
  
  #cbind(N, N.hg, L.hg, L.agr)
  dt.stan <- list(N=N, N_hg=N.hg, L_hg=L.hg, L_agr=L.agr)
  dt.stan <-lapply(dt.stan, function(x) if(length(x)==0){x <- 0}else{x <- x}  )
  
  return(dt.stan)
}

make_stan_data_auto <- function(dt){
  ## make data
  N=sum(dt$n) # total number of observation
  
  # number of HG languegaes
  N.hg <- ( dt %>% filter(Subsistence.AUTOTYP=='HG') %>% summarise(sum(n)) %>% c())
  N.hg <- N.hg$`sum(n)`
  
  # number of HG languages with labiodentals
  L.hg <- dt %>% filter(Subsistence.AUTOTYP=='HG', has.labiodentals==T) %>% c()
  L.hg <-L.hg$n
  
  # number of AGR languages with labiodentals
  L.agr <- dt %>% filter(Subsistence.AUTOTYP=='AGR', has.labiodentals==T) %>% c()
  L.agr <-L.agr$n
  
  #cbind(N, N.hg, L.hg, L.agr)
  dt.stan <- list(N=N, N_hg=N.hg, L_hg=L.hg, L_agr=L.agr)
  dt.stan <-lapply(dt.stan, function(x) if(length(x)==0){x <- 0}else{x <- x}  )
  
  return(dt.stan)
}