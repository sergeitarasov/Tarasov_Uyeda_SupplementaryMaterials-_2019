#==================
#
# Main plot for Causal Models
#
# 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
# Tarasov, Uyeda 2019
#==================


#==================
#
# Make dataset for plotting distrs
#
#==================
param.distr <-c()
for (i in 1:nrow(dt.freq)){
  block <- as.matrix(fit.tom.ind[[i]])
  block <- as_tibble(block)
  block <- add_column(block, Area=names(fit.tom.ind[i]) )
  block <-block[1:20000,]
  param.distr <-rbind(param.distr, block)
}
# remove Africa and NC Asia
names(fit.tom.ind)
param.distr <-param.distr %>% dplyr::filter(Area!='N-C Asia' & Area!='Africa' & Area!="Africa: Atlantic-Congo" & 
                                              Area!="Africa: non-Atlantic-Congo" & Area!="Whole World")
####
# dep model

block1 <- as.matrix(fit.tom.dep$`N-C Asia`)[1:20000,] %>% as_tibble() %>% add_column(Area='N-C Asia')
block2 <- as.matrix(fit.tom.dep$Africa)[1:20000,] %>% as_tibble() %>% add_column(Area='Africa')
As.Afr <- bind_rows(block1, block2)

param.distr <-bind_rows(As.Afr, param.distr)
param.distr %>% dplyr::filter(Area=='Papua')

####
# order Areas
ord <- param.distr %>% group_by(Area) %>% summarise(mean(theta_L1)) 
ord <- ord %>% arrange( `mean(theta_L1)`)
ord <- ord[[1]] #%>% rev


param.distr$Area <- as.factor(param.distr$Area)
levels(param.distr$Area )
#new.f <- factor(param.distr$Area,levels(param.distr$Area)[c(2, 4, 8, 7, 5, 6, 3, 9, 10, 1)])
new.f <- factor(param.distr$Area, ord)
levels(new.f  )
param.distr$Area <- new.f



#==================
#
# Make BF data
#
#==================

BF.melt <- BF[1:10,]

# order
BF.melt$Areas <- as.factor(BF.melt$Areas)
levels(param.distr$Area )
new.bf <- factor(BF.melt$Areas, ord)
levels(new.bf  )
BF.melt$Areas  <- new.bf

p.bf <- ggplot(BF.melt, aes(Areas)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=2.3, linetype="dashed", color = "black", size=0.3) +
  geom_hline(yintercept=-2.3, linetype="dashed", color = "black", size=0.3) +
  geom_vline(xintercept=c(1:10),  color = "black", size=0.1) +
  geom_point(aes(y=Tom),size=4, color="black", alpha=1, fill = "#E69F00", stroke = .1, shape = 21 ) +
  geom_point(aes(y=Auto),size=2,color="black", alpha=1, fill = "#56B4E9", stroke = .1, shape = 21 ) +
  theme(axis.text.x = element_text( size=12), axis.ticks.length=unit(.25, "cm"),  axis.ticks.y = element_blank(), axis.line.y =element_blank()  )+
  labs(x = NULL, y=NULL)+
  scale_x_discrete(labels = NULL)+
  coord_flip()+
  scale_y_continuous( limits=c(-6.2, 6.2), breaks=c(0, 1.16, 2.3, -1.16, -2.3, -4.6, 4.6), position = "left", labels = NULL)

#ggsave(paste0('BF1', '.eps'))

#==================
#
# HF gradient
#
#==================
library(ggridges)

Cont <- plyr::revalue(ord, c("W and SW Eurasia"="W Eurasia"))

p.hg <- 
  ggplot(
  param.distr, 
  aes(x = `theta_S`, y = `Area`)
  ) +
  geom_density_ridges_gradient(gradient_lwd=0.9, rel_min_height = 0.0,
                               aes(fill = ..x..), scale = 5, size = 0.1
  ) +
  scale_fill_gradientn(colours = terrain.colors(10) %>% rev)+
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), labels=Cont) +
  labs(title = '', y=NULL, x=NULL) +
  guides(fill=FALSE)+
  theme_ridges(font_size = 9, grid = F) + theme(axis.title.y = element_blank(), axis.text.x = element_text( size=6))

#==================
#
# LB gradient
#
#==================
library(RColorBrewer)

# set colors
cols1 <- brewer.pal(9,"Blues")

p.lb <- 
  ggplot(
  param.distr, 
  aes(x = `theta_L1`, y = `Area`)
) +
  geom_density_ridges_gradient(gradient_lwd=0.9, rel_min_height = 0.00,
                               aes(fill = ..x..), scale = 5, size = 0.1
  ) +
  scale_fill_gradientn(colours = cols1[2:6])+
  
  geom_density_ridges_gradient(gradient_lwd=0.9, rel_min_height = 0.00,
                               aes(x=`theta_L2`, fill = ..x..), scale =.4, size = 0.1
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  labs(title = '', y=NULL, x=NULL) +
  guides(fill=FALSE)+
  theme_ridges(font_size = 9, grid = F) + theme(axis.title.y = element_blank(), axis.text.x = element_text( size=6), plot.margin = margin(2.48, 2, 0.16, -.1, "cm"))

#==================
#
# Combine plots
#
#==================
library(gridExtra)

p.bf <- p.bf+theme(plot.margin = margin(2.7, 0, 0.1, 2, "cm"))
p.hg <- p.hg+theme(plot.margin = margin(0.1, 0.3, 0.3, .7, "cm"))
p.lb <- p.lb+theme(plot.margin = margin(1.63, 2, 0.3, -.1, "cm"))

grid.arrange(p.bf, p.hg, p.lb, nrow=1)

#library(gridExtra)
dev.off()
dev.new(width=160, height=100, unit="mm")
f1 <- grid.arrange(p.bf, p.hg, p.lb, nrow=1)
# ggsave('Fig1.png', f1, width = 160, height = 100, units = c('mm') )
# ggsave('Fig1.eps', f1, width = 160, height = 100, units = c('mm') )
