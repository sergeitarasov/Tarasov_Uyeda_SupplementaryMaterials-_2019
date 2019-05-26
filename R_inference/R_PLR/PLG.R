###########################
#
# Poisson Linear Regression (PLR): comparison of models with and without the Subsistence predictor
# Comment 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################
library("MASS")

#---- set directory
# set working directory to 'Tarasov_Uyeda_Supplementary' then
setwd("./R_inference/R_PLR")

# Read necessary objects and data ------------

# phoible.area (phoible_distance_dat.csv) is the same as the original phoible but includes distance from Africa
# this distance is obtained as described in PP_get_geo-dist.R
phoible.area<-read.csv("../R_PPS/output/phoible_distance_dat.csv", header=T, stringsAsFactors = F, sep=',')
phoible.area<-as_tibble(phoible.area)
phoible.area

#------ make datatset 
library("plyr")
phoible.area <-phoible.area %>% dplyr::mutate(Afr.HG= replace(rep(0,n()), which(Area=='Africa' & Subsistence.TOM=='HG'), 1 )  )
phoible.area <- phoible.area %>% dplyr::mutate(NC.HG= replace(rep(0,n()), which(Area=='N-C Asia' & Subsistence.TOM=='HG'), 1 )  )

# calculate the residual phonemes NonLb.NonFr.phonemes as (Nonlabiodental.phonemes-Nonlabiodental.Fricatives)
phoible.area <- phoible.area %>% dplyr::mutate(NonLb.NonFr.phonemes=Nonlabiodental.phonemes-Nonlabiodental.Fricatives)
phoible.area <- phoible.area %>% dplyr::mutate(Subsistence.TOM.num= mapvalues(Subsistence.TOM, c('AGR', 'HG'), c(0,1) ) %>% as.numeric() )



#----------------
#
# TOM dataset
#
#----------------
# filter data for PP simulations using TOM --------------
phoible.area.sel <- phoible.area %>% dplyr::select(Labiodentals,  Nonlabiodental.Fricatives, Area, Dist, 
                                                   Subsistence.TOM,  NonLb.NonFr.phonemes)

# remove 'Africa' and 'N-C Asia' since they support dependent model in Causal Models inference
phoible.area.sel <- phoible.area.sel %>% filter(Area!='Africa' & Area!='N-C Asia')
phoible.area.sel <-phoible.area.sel%>% na.omit()

#-------
# N0 Subsistence.TOM predictor
#-------

## Vectors of outcome and predictors
outcome    <- c("Labiodentals")

predictors <- c("Nonlabiodental.Fricatives", "NonLb.NonFr.phonemes", "Area", "Dist",
                'Nonlabiodental.Fricatives:Area', 'Nonlabiodental.Fricatives:NonLb.NonFr.phonemes', 'Nonlabiodental.Fricatives:Dist',
                'NonLb.NonFr.phonemes:Area', 'NonLb.NonFr.phonemes:Dist',
                'Dist:Area'
)


# List of models
list.models <- lapply(seq_along((predictors)), function(n) {
  
  rhs <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")
  
  paste(outcome, rhs, sep = "  ~  ")
})

## Vector of models
vector.of.models <- unlist(list.models)

## Fit models
fit.tb <- tibble(AIC=rep(NA, length(vector.of.models) ), Model=rep(NA, length(vector.of.models) ))
for (i in seq_along(vector.of.models))
  {
    formula <- as.formula(vector.of.models[i])
    fit <- glm(formula, data=phoible.area.sel, family="poisson")
    fit.tb[i,1] <- fit$aic
    fit.tb[i,2] <- vector.of.models[i]
  }

fit.Tom.noSub <- fit.tb %>% arrange(-desc(AIC))


#-------
# WITH Subsistence.TOM predictor
#-------

## Vectors of outcome and predictors
outcome    <- c("Labiodentals")

predictors <- c("Nonlabiodental.Fricatives", "NonLb.NonFr.phonemes", "Area", "Dist", 'Subsistence.TOM',
                'Nonlabiodental.Fricatives:Area', 'Nonlabiodental.Fricatives:NonLb.NonFr.phonemes', 'Nonlabiodental.Fricatives:Dist', 'Nonlabiodental.Fricatives:Subsistence.TOM',
                'NonLb.NonFr.phonemes:Area', 'NonLb.NonFr.phonemes:Dist', 'NonLb.NonFr.phonemes:Subsistence.TOM',
                'Dist:Area', 'Dist:Subsistence.TOM', 'Area:Subsistence.TOM'
)



# List of models
list.models <- lapply(seq_along((predictors)), function(n) {
  
  rhs <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")
  
  paste(outcome, rhs, sep = "  ~  ")
})

## Vector of models
vector.of.models <- unlist(list.models)

## Fit models
fit.tb <- tibble(AIC=rep(NA, length(vector.of.models) ), Model=rep(NA, length(vector.of.models) ))
for (i in seq_along(vector.of.models))
{
  formula <- as.formula(vector.of.models[i])
  fit <- glm(formula, data=phoible.area.sel, family="poisson")
  fit.tb[i,1] <- fit$aic
  fit.tb[i,2] <- vector.of.models[i]
}

fit.Tom.Sub <- fit.tb %>% arrange(-desc(AIC))

# save models in RDS
#saveRDS(fit.Tom.Sub, file = "output/PLR_all_models_TOM.rds")

#---- dAIC: Difference between best sub and NoSub models
(fit.Tom.noSub$AIC[1]-fit.Tom.Sub$AIC[1]) %>% round(., 2)
# 2.30

#-----

#--- Akaike weights
fit.Tom.Sub[2000,]
minAIC <- fit.Tom.Sub$AIC[1]
fit.Tom.Sub <-fit.Tom.Sub %>% mutate(dAIC=AIC-minAIC)
fit.Tom.Sub <-fit.Tom.Sub %>% mutate(AICw=exp(-0.5*dAIC)/sum(exp(-0.5*dAIC)) )
fit.Tom.Sub <-fit.Tom.Sub %>% mutate(AICw.cum=cumsum(AICw))
upp.lim <- which(fit.Tom.Sub$AICw.cum<.95)
upp.lim <-upp.lim[length(upp.lim)]
upp.lim # Numbwer of models within 0.95: 1892
fit.Tom.Sub[upp.lim,]

mm <- fit.Tom.Sub$Model[1]
grepl('Subsistence.TOM', mm)
sub.pr <- lapply(fit.Tom.Sub$Model[1:upp.lim], function(x) grepl('Subsistence.TOM', x)) %>% unlist

# subsitence present
which(sub.pr==T) %>% length() # Subsitence predictor is present in 1786 (out of 1892 models) 
(which(sub.pr==T) %>% length()/upp.lim) %>% round(., 3) # Subsitence predictor is present in 94%

(which(sub.pr==F) %>% length()) # Subsitence predictor is ABSENT in 106 models
(which(sub.pr==F) %>% length()/upp.lim) %>% round(., 3) # Subsitence predictor is ABSENT in 6%


# subsitence present without interactions
sub.prR <- lapply(fit.Tom.Sub$Model[1:1736], function(x) grepl('Subsistence.TOM:', x)) %>% unlist
rr <- which(sub.prR==T) %>% length()
sub.prL <- lapply(fit.Tom.Sub$Model[1:1736], function(x) grepl(':Subsistence.TOM', x)) %>% unlist
ll <- which(sub.prL==T) %>% length()
rr+ll # Subsitence  predictor is present as a predictor with interactions in 1523 models 
( (rr+ll)/upp.lim ) %>% round(., 2) # Subsitence.Tom predictor is present as a predictor with interactions in (80%) models

(upp.lim-((which(sub.pr==F) %>% length())) )- rr-ll# Subsitence.Tom predictor is present as a predictor without interactions in 263 (7%)
(((upp.lim-((which(sub.pr==F) %>% length())) )- rr-ll ) /upp.lim) %>% round(., 2) #Subsitence.Tom predictor is present as a predictor interactions in 14%

##########################################################################################################################################






#----------------
#
# AUTOTYP data
#
#----------------

phoible.area.sel <- phoible.area %>% dplyr::select(Labiodentals,  Nonlabiodental.Fricatives, Area, Dist, 
                                                   Subsistence.AUTOTYP, NonLb.NonFr.phonemes)
phoible.area.sel <- phoible.area.sel %>% filter(Area!='Africa' & Area!='N-C Asia')
phoible.area.sel <- phoible.area.sel %>%na.omit()

#-------
# N0 Subsistence.AUTO predictor
#-------

## Vectors of outcome and predictors
outcome    <- c("Labiodentals")

predictors <- c("Nonlabiodental.Fricatives", "NonLb.NonFr.phonemes", "Area", "Dist",
                'Nonlabiodental.Fricatives:Area', 'Nonlabiodental.Fricatives:NonLb.NonFr.phonemes', 'Nonlabiodental.Fricatives:Dist',
                'NonLb.NonFr.phonemes:Area', 'NonLb.NonFr.phonemes:Dist',
                'Dist:Area'
)


# List of models
list.models <- lapply(seq_along((predictors)), function(n) {
  
  rhs <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")
  
  paste(outcome, rhs, sep = "  ~  ")
})

## Vector of models
vector.of.models <- unlist(list.models)

## Fit models
fit.tb <- tibble(AIC=rep(NA, length(vector.of.models) ), Model=rep(NA, length(vector.of.models) ))
for (i in seq_along(vector.of.models))
{
  formula <- as.formula(vector.of.models[i])
  fit <- glm(formula, data=phoible.area.sel, family="poisson")
  fit.tb[i,1] <- fit$aic
  fit.tb[i,2] <- vector.of.models[i]
}

fit.AUTO.noSub <- fit.tb %>% arrange(-desc(AIC))


#-------
# WITH Subsistence.AUTO predictor
#-------

## Vectors of outcome and predictors
outcome    <- c("Labiodentals")


predictors <- c("Nonlabiodental.Fricatives", "NonLb.NonFr.phonemes", "Area", "Dist", 'Subsistence.AUTOTYP',
                'Nonlabiodental.Fricatives:Area', 'Nonlabiodental.Fricatives:NonLb.NonFr.phonemes', 'Nonlabiodental.Fricatives:Dist', 'Nonlabiodental.Fricatives:Subsistence.AUTOTYP',
                'NonLb.NonFr.phonemes:Area', 'NonLb.NonFr.phonemes:Dist', 'NonLb.NonFr.phonemes:Subsistence.AUTOTYP',
                'Dist:Area', 'Dist:Subsistence.AUTOTYP', 'Area:Subsistence.AUTOTYP'
)




# List of models
list.models <- lapply(seq_along((predictors)), function(n) {
  
  rhs <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")
  
  paste(outcome, rhs, sep = "  ~  ")
})

## Vector of models
vector.of.models <- unlist(list.models)

## Fit models
fit.tb <- tibble(AIC=rep(NA, length(vector.of.models) ), Model=rep(NA, length(vector.of.models) ))
for (i in seq_along(vector.of.models))
{
  formula <- as.formula(vector.of.models[i])
  fit <- glm(formula, data=phoible.area.sel, family="poisson")
  fit.tb[i,1] <- fit$aic
  fit.tb[i,2] <- vector.of.models[i]
}

fit.AUTO.Sub <- fit.tb %>% arrange(-desc(AIC))

# save models as RDS
#saveRDS(fit.AUTO.Sub, file = "output/PLR_all_models_AUTO.rds")


#--- Akaike weights
#fit.AUTO.Sub <- fit.AUTO.Sub %>% na.omit()
minAIC <- fit.AUTO.Sub$AIC[1]
fit.AUTO.Sub <-fit.AUTO.Sub %>% mutate(dAIC=AIC-minAIC)
fit.AUTO.Sub <-fit.AUTO.Sub %>% mutate(AICw=exp(-0.5*dAIC)/sum(exp(-0.5*dAIC)) )
fit.AUTO.Sub <-fit.AUTO.Sub %>% mutate(AICw.cum=cumsum(AICw))

upp.lim <-which(fit.AUTO.Sub$AICw.cum<.95)
upp.lim <-upp.lim[length(upp.lim)]
upp.lim # Numbwer of models within 0.95: 4549
sum(fit.AUTO.Sub$AICw[1:4549])

mm <- fit.AUTO.Sub$Model[1]
grepl('Subsistence.AUTOTYP', mm)
sub.pr <- lapply(fit.AUTO.Sub$Model[1:upp.lim], function(x) grepl('Subsistence.AUTOTYP', x)) %>% unlist
# models 1:6 do not contain Subsistence.AUTOTYP; 
# dAIC NoSub-Sub
fit.AUTO.Sub$AIC[1]-fit.AUTO.Sub$AIC[7]
#- 0.98 in favor of NoSub


# subsitence present
which(sub.pr==T) %>% length() # Subsitence predictor is present in 4127 (out of 4549 models) 
(which(sub.pr==T) %>% length()/upp.lim) %>% round(., 3) # Subsitence predictor is present in 0.907

(which(sub.pr==F) %>% length()) # Subsitence predictor is ABSENT in 422 models
(which(sub.pr==F) %>% length()/upp.lim) %>% round(., 3) # Subsitence predictor is ABSENT in 0.093


# subsitence present without interactions
sub.prR <- lapply(fit.AUTO.Sub$Model[1:4549], function(x) grepl('Subsistence.AUTOTYP:', x)) %>% unlist
rr <-which(sub.prR==T) %>% length()
sub.prL <- lapply(fit.AUTO.Sub$Model[1:4549], function(x) grepl(':Subsistence.AUTOTYP', x)) %>% unlist
ll <- which(sub.prL==T) %>% length()
rr+ll  # Subsitence.Tom predictor is present as a predictor with interactions in 3783 models (83%)
( (rr+ll)/upp.lim ) %>% round(., 2) # Subsitence.Tom predictor is present as a predictor with interactions in 0.83

(upp.lim-((which(sub.pr==F) %>% length())) )- rr-ll# Subsitence.Tom predictor is present as a predictor without interactions in 344
(((upp.lim-((which(sub.pr==F) %>% length())) )- rr-ll ) /upp.lim) %>% round(., 2) #Subsitence.Tom predictor is present as a predictor interactions in 8%

#################

