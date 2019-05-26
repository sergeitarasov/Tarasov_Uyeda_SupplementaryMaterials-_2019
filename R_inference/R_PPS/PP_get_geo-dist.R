###########################
#
# Estimating distance from Africa fro languages in phible dataset
#
# Comment on 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################

# set working directory to 'Tarasov_Uyeda_Supplementary' then
setwd("./R_inference/R_PP")

library(geosphere)

# get world map
map.world <- map_data("world")
# read data
phoible<-read.csv("../data/data_labiodentals_worldwide.csv",header=T, stringsAsFactors = F, sep=';')
phoible<-as_tibble(phoible)
phoible

#------------------------------------------------------------------
#
# Create set of geographic points that reflect human dispersal from Africa
#
# We follow the methodology of Atkinson (2011) but adjust the points
#  to the geographic Areas of the phoible dataset
#------------------------------------------------------------------

# Here is the set of the migrational geo-points
p.focal <-rbind(
  c(NA, NA),      # Origin
  c(32.078490, 30.027482), # from Africa
  c(116.803713,  -0.772777), # to Oceania
  c(-167.739947,  65.614918), # to N America
  c(-77.691820,   7.773037) # to S America
)

#------ find the best origin point
#
# Find the best origin point in Africa 
# that best fit the decrease of phonemes using linear regression
#------

phoible.area1 <- phoible %>% filter(Area=="Africa" & Longitude!=is.na(Longitude))
nonLab <- list()
for (i in 1:nrow( phoible.area1)){
  point <- c(phoible.area1$Longitude[i],  phoible.area1$Latitude[i])
  dist <- apply(phoible.area1, 1, function(x) distGeo(point, as.numeric(c(x[13], x[12]) ) )/1000  )
  phoible.area1 <- bind_cols(phoible.area1, Dist=dist)
  nonLab[[i]] <- glm(Nonlabiodental.phonemes ~ Dist, data=phoible.area1, family="poisson") 
  phoible.area1 <- phoible.area1 %>% select(-Dist)
}

aic <- lapply(nonLab, function(x) AIC(x)) %>% unlist
best.aic <- which(aic==min(aic) )
nonLab[best.aic]
phoible.area1[best.aic,]
# the best point is the language kgal1244
best.coord <- phoible.area1[best.aic,c('Longitude', 'Latitude')]
# Coords: 22.7    -25.9

# plot position of kgal1244
ggplot(phoible.area1[best.aic, ]) + 
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group), fill='white', color='white', size=.15)+
  geom_point( aes(x = Longitude, y = Latitude), alpha=1, size=2 )+
  coord_equal()

#------

# Finilize the the set of the migrational geo-points
p.focal <-rbind(
  best.coord,      # Origin
  c(32.078490, 30.027482), # from Africa
  c(116.803713,  -0.772777), # to Oceania
  c(-167.739947,  65.614918), # to N America
  c(-77.691820,   7.773037) # to S America
)

#colnames(p.focal) <- c('Longitude', 'Latitude')
#p.focal <-as_tibble(p.focal)

perArea <- list(
  Origin=c('Africa'),
  Afr=c("W and SW Eurasia", "S/SE Asia", "N-C Asia"),
  Oceania=c("Papua", "Pacific", "Australia"),
  NAmer=c("N America", "C America"),
  SAmer=c("S America")
)

# Calculate the distance (km) to language origin in Africa using point from p.focal
  Origin=c(0)
  Afr=distGeo(p.focal[1,], c(32.078490, 30.027482))/1000
  Oceania=Afr+distGeo(c(32.078490, 30.027482), c(116.803713,  -0.772777) )/1000
  NAmer=Afr+distGeo(c(32.078490, 30.027482),c(-167.739947,  65.614918))/1000
  SAmer=NAmer+distGeo(c(-167.739947,  65.614918), c(-77.691820,   7.773037))/1000

d2org <- c(Origin, Afr, Oceania, NAmer, SAmer)

# Plot all migrational points
ggplot(p.focal) + 
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group), fill='white', color='white', size=.15)+
  geom_point( aes(x = Longitude, y = Latitude), alpha=1, size=2 )+
  coord_equal() +
  theme(
    panel.background = element_rect(fill = 'grey', linetype = 0)
    ,panel.grid = element_blank()
    ,plot.title = element_text(size = 30)
    ,plot.subtitle = element_text(size = 10)
    ,axis.text = element_blank()
    ,axis.title = element_blank()
    ,axis.ticks = element_blank()
    ,axis.line = element_blank()
    ,axis.line.x = element_blank()
)

# Plot language geographic arear
ggplot(phoible) + 
  #ggplot(phoible %>% filter(Area=='W and SW Eurasia')) + 
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group), fill='white', color='white', size=.15)+
  geom_point( aes(x = Longitude, y = Latitude, color=Area), alpha=1, size=1 )+
  geom_point(data =p.focal, aes(x = Longitude, y = Latitude), alpha=1, size=2 )+
  coord_equal() +
  theme(
    panel.background = element_rect(fill = 'grey', linetype = 0)
    ,panel.grid = element_blank()
    ,plot.title = element_text(size = 30)
    ,plot.subtitle = element_text(size = 10)
    ,axis.text = element_blank()
    ,axis.title = element_blank()
    ,axis.ticks = element_blank()
    ,axis.line = element_blank()
    ,axis.line.x = element_blank()
  )


#----------------------------------------------------
#
# Estimate distance from the origin to all languages
#
#----------------------------------------------------

phoible.area1 <- phoible %>% filter(Area=="Africa")
dist <- apply(phoible.area1, 1, function(x) distGeo(p.focal[1,], as.numeric(c(x[13], x[12]) ) )/1000  )
phoible.area1 <- bind_cols(phoible.area1, Dist=dist)

phoible.area2 <- phoible %>% filter(Area=="W and SW Eurasia" | Area=="S/SE Asia" | Area=="N-C Asia" )
dist <- d2org[2] + apply(phoible.area2, 1, function(x) distGeo(p.focal[2,], as.numeric(c(x[13], x[12]) ) )/1000  )
phoible.area2 <- bind_cols(phoible.area2, Dist=dist)
plot(phoible.area2$Dist, phoible.area2$Nonlabiodental.phonemes)

phoible.area3 <- phoible %>% filter(Area=="Papua" | Area=="Pacific" | Area=="Australia" )
dist <- d2org[3] + apply(phoible.area3, 1, function(x) distGeo(p.focal[3,], as.numeric(c(x[13], x[12]) ) )/1000  )
phoible.area3 <- bind_cols(phoible.area3, Dist=dist)
plot(phoible.area3$Dist, phoible.area3$Nonlabiodental.phonemes)

phoible.area4 <- phoible %>% filter(Area=="N America" | Area=="C America")
dist <- d2org[4] + apply(phoible.area4, 1, function(x) distGeo(p.focal[4,], as.numeric(c(x[13], x[12]) ) )/1000  )
phoible.area4 <- bind_cols(phoible.area4, Dist=dist)
plot(phoible.area4$Dist, phoible.area4$Nonlabiodental.phonemes)

phoible.area5 <- phoible %>% filter(Area=="S America")
dist <- d2org[5] + apply(phoible.area5, 1, function(x) distGeo(p.focal[5,], as.numeric(c(x[13], x[12]) ) )/1000  )
phoible.area5 <- bind_cols(phoible.area5, Dist=dist)
plot(phoible.area5$Dist, phoible.area5$Nonlabiodental.phonemes)

phoible.area <- bind_rows(phoible.area1, phoible.area2, phoible.area3, phoible.area4, phoible.area5)

#------ save file for further use
#saveRDS(phoible.area, file = 'output/phoible.area.rds')
#write.csv(phoible.area, file = 'output/phoible_distance_dat.csv')
