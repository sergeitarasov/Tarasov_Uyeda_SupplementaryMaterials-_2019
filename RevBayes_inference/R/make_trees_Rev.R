
# Read in data
labiodental.data <- read.csv('ie_labiodental_data_d.csv', stringsAsFactors=F)

# get list of prsent languages
glot <- paste0('<', labiodental.data$glottocode, '>')
Langs <- paste0(labiodental.data$language, glot)

#----------- create tre files by sampling 1000 trees

# ----- C dataset
ie.c.trees <- read.nexus('indoeuropean-C-20k-matched.nex')
ie.c.trees <- lapply(ie.c.trees, function(t) {
  t <- drop.tip(t, setdiff(t$tip.label, Langs))
  return(t)
})
class(ie.c.trees) <- 'multiPhylo'	
attributes(ie.c.trees)$TipLabel <- ie.c.trees[[1]]$tip.label

# sample 1000 trees
sample <- sample(c(1:length(ie.c.trees)), 1000)
tree.sampled <- ie.c.trees[sample]
write.tree(tree.sampled, file = 'indoeuropean-C.tre')



# ----- B dataset
ie.c.trees <- read.nexus('indoeuropean-B-10k-matched.nex')
ie.c.trees <- lapply(ie.c.trees, function(t) {
  t <- drop.tip(t, setdiff(t$tip.label, Langs))
  return(t)
})
class(ie.c.trees) <- 'multiPhylo'	
attributes(ie.c.trees)$TipLabel <- ie.c.trees[[1]]$tip.label

# sample 1000 trees
sample <- sample(c(1:length(ie.c.trees)), 1000)
tree.sampled <- ie.c.trees[sample]
write.tree(tree.sampled, file = 'indoeuropean-B.tre')

