###########################
#
# Fitting dependent and independent binomial causal models for the language dataset from
# 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################

# installing Stan correctly

remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
pkgbuild::has_build_tools(debug = TRUE)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native",
    if( grepl("^darwin", R.version$os)) "CXX14FLAGS += -arch x86_64 -ftemplate-depth-256" else 
      if (.Platform$OS.type == "windows") "CXX11FLAGS=-O3 -march=native -mtune=native" else
        "CXX14FLAGS += -fPIC",
    file = M, sep = "\n", append = TRUE)

library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)

#########
# libraries

library("rstan")
library("dplyr")
library("tibble")
library("viridis")
library('bridgesampling')

##########################
#
# Read in Phoible dataset
#

# set working directory to 'Tarasov_Uyeda_Supplementary' then
setwd("./R_inference/data")

phoible<-read.csv("data_labiodentals_worldwide.csv",header=T, stringsAsFactors = F, sep=';')
phoible<-as_tibble(phoible)
phoible

######
#
# source useful functions for making Stan data
setwd("../R_CM")
source(file='CM_func_make_stan_data.R')

###########
#
# Make datasets for the analyses: TOM
#
###########

# get unique areas
Areas <- phoible %>% select( Area ) %>%unique() 


dt.freq <- c()
for (i in 1:nrow(Areas)){
  dt <- phoible %>% filter( Area==Areas[[i,1]] ) %>% group_by(Subsistence.TOM) %>% count(has.labiodentals)
  dt <- make_stan_data(dt) %>% as_tibble()
  dt.freq <- rbind(dt.freq, dt)
}
dt.freq <- add_column(dt.freq, Areas=Areas[[1,]], .before=1)


# all world
dt <- phoible %>% group_by(Subsistence.TOM) %>% count(has.labiodentals)
dt <- make_stan_data(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Whole World', .before=1)
dt.freq <- bind_rows(dt.freq, dt ) 

# Africa: Africaca Atlantic-Congo vs. Africa the rest
# Most languages (422 langs) in Africa belong to Atlantic-Congo family, let's see the summary
phoible %>% filter( Area=='Africa' ) %>% group_by(Family) %>% tally(sort = T)

# separate Atl-Congo
dt <- phoible %>% filter( Area=='Africa', Family=='Atlantic-Congo [atla1278]' ) %>%
  group_by(Subsistence.TOM) %>% count(has.labiodentals)
dt <- make_stan_data(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Africa: Atlantic-Congo', .before=1)
dt.freq <- bind_rows(dt.freq, dt ) 

# separate non-Atl-Congo
dt <- phoible %>% filter( Area=='Africa', Family!='Atlantic-Congo [atla1278]' ) %>%
  group_by(Subsistence.TOM) %>% count(has.labiodentals)
dt <- make_stan_data(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Africa: non-Atlantic-Congo', .before=1)
dt.freq <- bind_rows(dt.freq, dt )

#-- 
# For the interpretation of columns in dt.freq see file CM_stan_models.R
#--


###########
#
# Make datasets: AUTO
#
###########
Areas <- phoible %>% select( Area ) %>%unique() 


dt.freq.auto <- c()
for (i in 1:nrow(Areas)){
  dt <- phoible %>% filter( Area==Areas[[i,1]] ) %>% group_by(Subsistence.AUTOTYP) %>% count(has.labiodentals)
  dt <-dt %>% filter( Subsistence.AUTOTYP=='AGR' | Subsistence.AUTOTYP=='HG')
  dt <- make_stan_data_auto(dt) %>% as_tibble()
  dt.freq.auto <- rbind(dt.freq.auto, dt)
}

dt.freq.auto <- add_column(dt.freq.auto, Areas=Areas[[1,]], .before=1)

# all world
dt <- phoible %>% group_by(Subsistence.AUTOTYP) %>% count(has.labiodentals)
dt <-dt %>% filter( Subsistence.AUTOTYP=='AGR' | Subsistence.AUTOTYP=='HG')
dt <- make_stan_data_auto(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Whole World', .before=1)
dt.freq.auto <- bind_rows(dt.freq.auto, dt ) 

# Africa: Africaca Atlantic-Congo vs. Africa the rest
# separate Atl-Congo
dt <- phoible %>% filter( Area=='Africa', Family=='Atlantic-Congo [atla1278]' ) %>%
  group_by(Subsistence.AUTOTYP) %>% count(has.labiodentals)
dt <-dt %>% filter( Subsistence.AUTOTYP=='AGR' | Subsistence.AUTOTYP=='HG')
dt <- make_stan_data_auto(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Africa: Atlantic-Congo', .before=1)
dt.freq.auto <- bind_rows(dt.freq.auto, dt )

# separate non-Atl-Congo
dt <- phoible %>% filter( Area=='Africa', Family!='Atlantic-Congo [atla1278]' ) %>%
  group_by(Subsistence.AUTOTYP) %>% count(has.labiodentals)
dt <-dt %>% filter( Subsistence.AUTOTYP=='AGR' | Subsistence.AUTOTYP=='HG')
dt <- make_stan_data_auto(dt) %>% as_tibble()
dt <- dt %>% add_column(Areas='Africa: non-Atlantic-Congo', .before=1)
dt.freq.auto <- bind_rows(dt.freq.auto, dt )
##########





########################################################################
############
# BI inference: initialize stan models
#############
source(file='CM_stan_models.R')
stan.M_G1_dep_Log <- stan_model(model_code = M_G1_dep_Log)
stan.M_G1_ind_Log <- stan_model(model_code = M_G1_ind_Log)
stan.M_G1_dep <- stan_model(model_code = M_G1_dep)
stan.M_G1_ind <- stan_model(model_code = M_G1_ind)

############
# BI inference  TOM
#############

# Independemt model
fit.tom.ind <- vector(mode = "list", length = 13)
fit.tom.ind.ML <- vector(mode = "list", length = 13)
names(fit.tom.ind) <- dt.freq[[1,]]
names(fit.tom.ind.ML) <- dt.freq[[1,]]

for (i in 1:nrow(dt.freq)){
  dt.stan <- dt.freq[i,]
  dt.stan <-list(N=dt.stan[[1,2]], N_hg=dt.stan[[1,3]], L_hg=dt.stan[[1,4]], L_agr=dt.stan[[1,5]] )
  
  # Fit indep
  fit.tom.ind[[i]]<- sampling(stan.M_G1_ind_Log, data = dt.stan,
                      iter = 50000, warmup = 5000, chains = 3, cores = 1)
  # calculation ML
  fit.tom.ind.ML[[i]]<- bridge_sampler(fit.tom.ind[[i]], silent = FALSE)
}



# Dependemt model
fit.tom.dep <- vector(mode = "list", length = 13)
fit.tom.dep.ML <- vector(mode = "list", length = 13)
names(fit.tom.dep) <- dt.freq[[1,]]
names(fit.tom.dep.ML) <- dt.freq[[1,]]

for (i in 1:nrow(dt.freq)){
  dt.stan <- dt.freq[i,]
  dt.stan <-list(N=dt.stan[[1,2]], N_hg=dt.stan[[1,3]], L_hg=dt.stan[[1,4]], L_agr=dt.stan[[1,5]] )
  
  # Fitdep
  fit.tom.dep[[i]]<- sampling(stan.M_G1_dep_Log, data = dt.stan,
                              iter = 50000, warmup = 5000, chains = 3, cores = 1)
  # calculation ML
  fit.tom.dep.ML[[i]]<- bridge_sampler(fit.tom.dep[[i]], silent = FALSE)
}


ML.tom.ind <- lapply(fit.tom.ind.ML, function(x) x$logml) %>% unlist
ML.tom.dep <-lapply(fit.tom.dep.ML, function(x) x$logml) %>% unlist

round(ML.tom.ind-ML.tom.dep, 3)


############
# BI inference AUTOTYP
#############

# Independemt model
fit.auto.ind <- vector(mode = "list", length = 13)
fit.auto.ind.ML <- vector(mode = "list", length = 13)
names(fit.auto.ind) <- dt.freq.auto[[1,]]
names(fit.auto.ind.ML) <- dt.freq.auto[[1,]]

for (i in 1:nrow(dt.freq.auto)){
  dt.stan <- dt.freq.auto[i,]
  dt.stan <-list(N=dt.stan[[1,2]], N_hg=dt.stan[[1,3]], L_hg=dt.stan[[1,4]], L_agr=dt.stan[[1,5]] )
  
  # Fit indep
  fit.auto.ind[[i]]<- sampling(stan.M_G1_ind_Log, data = dt.stan,
                              iter = 50000, warmup = 5000, chains = 3, cores = 1)
  # calculation ML
  fit.auto.ind.ML[[i]]<- bridge_sampler(fit.auto.ind[[i]], silent = FALSE)
}



# Dependemt model
fit.auto.dep <- vector(mode = "list", length = 13)
fit.auto.dep.ML <- vector(mode = "list", length = 13)
names(fit.auto.dep) <- dt.freq.auto[[1,]]
names(fit.auto.dep.ML) <- dt.freq.auto[[1,]]

for (i in 1:nrow(dt.freq.auto)){
  dt.stan <- dt.freq.auto[i,]
  dt.stan <-list(N=dt.stan[[1,2]], N_hg=dt.stan[[1,3]], L_hg=dt.stan[[1,4]], L_agr=dt.stan[[1,5]] )
  
  # Fitdep
  fit.auto.dep[[i]]<- sampling(stan.M_G1_dep_Log, data = dt.stan,
                              iter = 50000, warmup = 5000, chains = 3, cores = 1)
  # calculation ML
  fit.auto.dep.ML[[i]]<- bridge_sampler(fit.auto.dep[[i]], silent = FALSE)
}


ML.auto.ind <- lapply(fit.auto.ind.ML, function(x) x$logml) %>% unlist
ML.auto.dep <-lapply(fit.auto.dep.ML, function(x) x$logml) %>% unlist

round(ML.auto.ind-ML.auto.dep, 3)


## BFs
cbind(ML.tom.ind-ML.tom.dep, ML.auto.ind-ML.auto.dep)
BF <- cbind(ML.tom.ind-ML.tom.dep, ML.auto.ind-ML.auto.dep) #%>% round(3)
colnames(BF) <- c('Tom', 'Auto')
BF <- as_tibble(BF)
BF <- add_column(BF, Areas=names(ML.tom.ind), .before=1)

#print(xtable(newobject2, type = "latex"), file = "filename2.tex")
# Latex table
library(xtable)
xtable(BF , type = "latex", digits = 3, caption ='BF')

