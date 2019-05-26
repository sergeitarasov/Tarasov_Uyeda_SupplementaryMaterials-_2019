#' ---
#' title: "The Evolution of Labiodental Articulation in Indo-European"
#' author: "Balthasar Bickel and Paul Widmer"
#' date: "July 2017" 
#' header-includes:
#'  - \usepackage[Symbol]{upgreek}
#'  - \usepackage{tocloft}
#'  - \settowidth{\cftsubsubsecnumwidth}{S4.5.3.10.}
#'  - \settowidth{\cftsubsecnumwidth}{S4.5.3.}
#'  - \settowidth{\cftsecnumwidth}{S4.5i.}
#'  - \renewcommand{\thefigure}{S4.\arabic{figure}}
#'  - \renewcommand{\thesection}{S4.\arabic{section}}
#' output:
#'  pdf_document:
#'   fig_crop: true
#'   fig_caption: true
#'   latex_engine: xelatex
#'   keep_tex: true
#'   toc: true
#'   toc_depth: 4
#'   number_section: true
#'   pandoc_args: [
#'   "--bibliography=/Users/bbickel/SWITCHdrive/Bibliographies/repo/bbbib.bib",
#'   "--csl=/Users/bbickel/SWITCHdrive/Bibliographies/repo/science.csl",
#'   ]
#' ---

#+ setup, include=F
library(parallel)
library(testthat)
library(devtools)
install_github("balthasarbickel/btw") # forked from "rgriff23/btw" with some additions
library(btw)
.BayesTraitsPath <- "/usr/local/bin/BayesTraitsV2"
library(ape)
install_github("liamrevell/phytools") # get the latest update, with a fix to simmap
library(phytools)
library(ggtree)
library(ggridges)
source_url('https://raw.githubusercontent.com/balthasarbickel/patches/f94801a3e9f1288f1117cd82ae3b21d67b841922/densityMap.R') # small patch removing color bar subtitle. Might not work with future releases of phytools!
source_url('https://raw.githubusercontent.com/balthasarbickel/patches/b25380247ba297a7c422e970e897ac689693170d/phenogram.R') # small patch of phytools' phenogram for more customized plotting
library(RColorBrewer)
library(phangorn)
library(ggplot2)
library(gridExtra)
library(ggraph)
library(igraph)
library(dplyr)
library(tidyr)
library(gdata)
library(knitr)
knit_hooks$set(crop = hook_pdfcrop, pars = function(before, options, envir) {if(before) {par(family=my.font)} else NULL}) # trick for keeping par across chunks; define my.font below!
opts_chunk$set(fig.path='figures/',
	       	dev = 'cairo_pdf', dev.args=list(bg='transparent'), # or quartz_pdf (for lattice)
                fig.height = 7,
                fig.width = 14,
                message = F,
                warning = F,
	        autodep=T,
	      	cache.comments=F,
	      	crop=T,
		pars=T
	      	)
# graphics setup:
my.font = 'Helvetica'   
# ggplot   
theme_set(theme_bw(base_size = 24) +
           theme(text = element_text(family = my.font),
                 plot.background = element_rect(fill = "transparent", colour = NA)
           ))

theme_phylo <- function(base_size=11, time=c('grid','timeline', 'none'), year.size=12, legend=NULL) {
	if(!is.null(legend)) {lt <- element_text(size=base_size)} else {lt <- element_blank()}
	switch(time,
		grid = {
		theme(axis.text.x=element_text(size=year.size), 
		  panel.grid.major.x=element_line(size=.2, color='grey', linetype='dashed'),
		  legend.position=legend, legend.title=lt
		  ) 
	  	},
		timeline= { 
		theme(axis.text.x=element_text(size=year.size),
				axis.ticks.length=unit(4,"mm"), axis.ticks=element_line(color='grey'),
		  legend.position=legend, legend.title=lt) 
	  	},
		none= {
		theme(axis.text.x=element_blank(), legend.position=legend, legend.title=lt)	
		}
		)
	}
reverse.time <- function(p) { p$data$x <- p$data$x - max(p$data$x); return(p) }

options(width=180)

#'   
#' 
#' 
#' \clearpage
#' 
#' Introduction
#' ===========
#' 
#' 
#' Indo-European comparative linguistics has identified what is called *correspondence sets*: words with related meanings show systematic sound correspondences. For example, there is a large set of related words that begin with bilabial *w* in some languages, e.g. in English (*will*, *word*, *wind* etc.). These word-initial sounds systematically correspond to a bilabially articulated sound in some related languages (e.g. Latin *uolō*, *uerbum*, *uentus* etc.) and to a labiodentally articulated sound in others (e.g. German *wollen*, *Wort*, *Wind* etc.). This establishes a set of correspondences that contains the articulations $\{w, v\}$ (among others). Each set reflects descent (with or without modification) from an ancestral proto-sound, but the reconstruction of the phonetic quality of this sound is an unresolved problem. Despite this problem, correspondence sets are traditionally labeled by a symbol that makes a vague suggestion on the phonetics, e.g. the $\{w, v\}$ set is conventionally labelled "$\text{*}w$" (with a distinctive asterisk), suggesting a bilabial realization. It is clear for experts, however, that these labels are mere suggestions and do not reflect an actual reconstruction of phonetics. To avoid confusion, we refer to the label of a correspondence set as a 'correspondeme'. Unlike a phoneme, a correspondeme is not associated with a well-defined range of phonetic values. A correspondeme only represents a distinct (meaning-differentiating) proto-sound from which all modern variants in a set derive, whatever that sound might have been phonetically. This observation challenges any inference from correspondeme labels to the presence or absence of labiodental places of articulation (an inference that was explicitly made by Hockett [@Hockett1985Distinguished]).
#'
#' Traditional attempts at reconstructing actual phonetic values of correspondemes  are chiefly based on qualitative considerations of maximum parsimony in the phylogeny, coupled with a tacit assumption that the values closely reflect what is reported in the earliest phonetic descriptions of well-investigated western languages (e.g. bilabial realization in Latin and Ancient Greek for decendants of the $\{w, v\}$ set). This approach incurs two problems. First, maximum-parsimony methods are problematic because sound changes can be completely reversed within relatively short time spans. For example the Germanic nasal $m$ became $f$ before $n$ in North Germanic (Old Norse, early Old Swedish *nafn* 'name') and reverted to $m$ in late Old Swedish (*namn*). Second, arguments based on earliest descriptions or even attestations are problematic because all current estimates show that the time that elapsed between Proto-Indo-European and the earliest attested languages is at least as long as the time between these languages and now. This is true even for relative young age estimates [@Changetal2015Ancestry-constrained].
#' 
#' Here we take a fresh approach to the problem and model sound change as Continuous-Time Markov Chains (CTMC), allowing for reversals and different transition rates between states (Section 3). We fit CTMC models using the Markov Chain Monte Carlo Sampling (MCMC) implemented in `BayesTraits` [@Pageletal2014BayesTraits]. We then take the best fitting models to estimate ancestral values (Section S4.4) and to reconstruct values for each time interval of Indo-European (Section S4.5). For this we apply Stochastic Character Mapping [@Nielsen2002Mapping;@Huelsenbecketal2003Stochastic;@Revell2012phytools], a method that has proven valid for linguistic reconstruction elsewhere [@Widmeretal2017NP-recursion].
#' 
#' 
#' 
#' 
#' Data and Phylogenies
#' ====================
#' 
#' Etymological research in Indo-European linguistics has established ten major correspondemes which have a labiodential reflex in at least one daughter language: $\text{*}b^h$, $\text{*}p$, $\text{*}w$, $\text{*}m$, $\text{*}g^w$, $\text{*}g^{wh}$, $\text{*}k^w$, $\text{*}d^h$, $\text{*}k^j$, and $\text{*}b$. The phonetic values of correspondemes are justified in detail in Section S6: 
#' 
#' 
#+ data, include=T, cache=T
labiodental.data <- read.csv('../data/ie_labiodental_data_d.csv', stringsAsFactors=F)
	
#+ extant, include=F, cache=F
correspondemes <- c('bh','p','w','m','gw','gwh','kw','dh','ky','b')
extant.languages <- apply(labiodental.data[labiodental.data$status %in% 'alive', correspondemes], 
                          2, function(x) { m <-  ifelse(x %in% '-', NA, x)})
#' Extant languages have `r round(mean(apply(extant.languages, 1, function(x) any(x %in% 'Y', na.rm=T)))*100,0)`% labiodentals. Figures S5.1 and S5.2 show their distribution in two different phylogenies of Indo-European (using the `ggtree`[@Yuetal2017ggtree]). We use phylogenies from Chang and colleagues[@Changetal2015Ancestry-constrained] (identified by the suffix `.c` here) and from Bouckaert and colleagues[@Bouckaertetal2012Mapping] (suffix `.b`):
#' 
#+ trees, cache=T
ie.c.sum.tree <- read.nexus('../data/indoeuropean-C-sum-matched.nex')
ie.c.sum.tree$tip.label <- gsub('(.*)(<.*>)','\\1', ie.c.sum.tree$tip.label)
ie.c.sum.tree <- drop.tip(ie.c.sum.tree, 
						setdiff(ie.c.sum.tree$tip.label, labiodental.data$language))

ie.b.sum.tree <- read.nexus('../data/indoeuropean-B-sum-matched.nex')
ie.b.sum.tree$tip.label <- gsub('(.*)(<.*>)','\\1', ie.b.sum.tree$tip.label)
ie.b.sum.tree <- drop.tip(ie.b.sum.tree, 
						setdiff(ie.b.sum.tree$tip.label, labiodental.data$language))

ie.c.trees <- read.nexus('../data/indoeuropean-C-20k-matched.nex')
ie.c.trees <- lapply(ie.c.trees, function(t) {
	t$tip.label <- gsub('(.*)(<.*>)','\\1', t$tip.label)
	t <- drop.tip(t, setdiff(t$tip.label, labiodental.data$language))
	return(t)
	})
	class(ie.c.trees) <- 'multiPhylo'	
	attributes(ie.c.trees)$TipLabel <- ie.c.trees[[1]]$tip.label

ie.b.trees <- read.nexus('../data/indoeuropean-B-10k-matched.nex')
ie.b.trees <- lapply(ie.b.trees, function(t) {
	t$tip.label <- gsub('(.*)(<.*>)','\\1', t$tip.label)
	t <- drop.tip(t, setdiff(t$tip.label, labiodental.data$language))
	return(t)
	})
	class(ie.b.trees) <- 'multiPhylo'	
	attributes(ie.b.trees)$TipLabel <- ie.b.trees[[1]]$tip.label
	
# /*
# setdiff(labiodental.data$language,ie.c.sum.tree$tip.label)
# setdiff(ie.c.sum.tree$tip.label,labiodental.data$language)
# setdiff(labiodental.data$language,ie.c.trees[[138]]$tip.label)
# setdiff(ie.c.trees[[138]]$tip.label,labiodental.data$language)
# setdiff(labiodental.data$language,ie.b.sum.tree$tip.label)
# setdiff(ie.b.sum.tree$tip.label,labiodental.data$language)
# setdiff(labiodental.data$language,ie.b.trees[[138]]$tip.label)
# setdiff(ie.b.trees[[138]]$tip.label,labiodental.data$language)
# */
#'  
#' 
#' 
#+ treedata-c,fig.height=14.5, fig.cap='Labiodental reflexes in each correspondeme that contains at least one labiodental reflex. Red indicates presence, blue absence a labiodental reflex. White indicates unknown or uncertain data. Maximum clade consensus tree from Chang et al.', cache=T, echo=F

ie.c.sum.tree.print <- ie.c.sum.tree
ie.c.sum.tree.print$tip.label <- gsub('_',' ',ie.c.sum.tree.print$tip.label)

ie.c.sum.tree.plot <- reverse.time(ggtree(ie.c.sum.tree.print, ladderize=T, right=T)) +
      geom_tiplab(size=4.5, align=T, linesize = .2)
	
labiodental.data.print <- labiodental.data	
rownames(labiodental.data.print) <- gsub('_',' ',labiodental.data.print$language)

color.scheme <-        c('white','blue', 'red')
names(color.scheme) <- c('-',    'N',     'Y')

colnames(labiodental.data.print)[4:13] <- c("bʰ","p","w","m","gʷ","gʷʰ","kʷ","dʰ","kʲ","b")
gheatmap(ie.c.sum.tree.plot, labiodental.data.print[,4:13],
	colnames_position='top', width=.4, offset=1250, color='black') +
		scale_fill_manual(name="", values=color.scheme) + 
		scale_x_continuous(breaks=c(-6000,-4000,-2000,0)) +
		scale_y_continuous(expand = c(-0.01, 1)) +
		theme_tree2(axis.text.x=element_text(size=12)) +
		theme(legend.position='none', axis.ticks.length=unit(4,"mm"), 
	 	 							  axis.ticks=element_line(color='grey'))	 
									  
# /*
# load them_phylo() above
# c.tree.plot <- rotate(ggtree(ie.c.sum.tree, ladderize=T, right=T),77) %>%
# rotate(87) %>% rotate(78) %>% rotate(82)
# c.treedata.plot <- c.tree.plot %<+% labiodental.data.print
# reverse.time(c.treedata.plot) +
# geom_tiplab(aes(label=label), offset=100) +
# 	# geom_label2(aes(x=x,y=y, label=node, subset=!isTip), size=2) +
# 	geom_tippoint(aes(fill=w), size=4, shape=22) +
# 		scale_fill_manual(name= '*w ', values=color.scheme,
# 						  labels=c('unknown', '-labiodental', '+labiodental')) +
# 		xlim(-6800,1200) +
# 			theme_phylo(time='none', legend='top') +
# 			theme(legend.text=element_text(size=20), legend.title=element_text(size=20))

c.treedata.plot <- ggtree(ie.c.sum.tree, ladderize=T, right=T) %<+% labiodental.data.print
c.treelist <- lapply(colnames(labiodental.data.print)[4:13], function(t) { reverse.time(c.treedata.plot) +
geom_tiplab(aes(label=label), offset=100, size=1.8) +
    geom_tippoint(aes_string(fill=t), size=2, shape=22) +
        scale_fill_manual(values=color.scheme) +
        xlim(-6800,1400) +
        annotate('text', x=-6500,y=26,label=paste("*",t)) + 
            theme_phylo(time='grid')
            })
do.call("grid.arrange", c(c.treelist, ncol=5))

# */
#' 
#' 
#' \clearpage
#' 
#+ treedata-b,fig.height=14.5, fig.cap='Labiodental reflexes in each correspondeme that contains at least one labiodental reflex. Red indicates presence, blue absence a labiodental reflex. White indicates unknown or uncertain data. Maximum clade consensus tree from Bouckaert et al', cache=T, echo=F

ie.b.sum.tree.print <- ie.b.sum.tree
ie.b.sum.tree.print$tip.label <- gsub('_',' ',ie.b.sum.tree.print$tip.label)

ie.b.sum.tree.plot <- reverse.time(ggtree(ie.b.sum.tree.print, ladderize=T, right=T)) +
      geom_tiplab(size=4.5, align=T, linesize = .2)
	
gheatmap(ie.b.sum.tree.plot, labiodental.data.print[,4:13],
	colnames_position='top', width=.4, offset=1500, color='black') +
		scale_fill_manual(name="", values=color.scheme) + 
		scale_x_continuous(breaks=c(-8000, -6000,-4000,-2000,0)) +
		scale_y_continuous(expand = c(-0.01, 1)) +
		theme_tree2(axis.text.x=element_text(size=12)) +
		theme(legend.position='none', axis.ticks.length=unit(4,"mm"), 
	 	 							  axis.ticks=element_line(color='grey'))	 

# /*
# for print (in main paper):
# quartz(height=17, width=11)
# offset=2000
b.treedata.plot <- ggtree(ie.b.sum.tree, ladderize=T, right=T) %<+% labiodental.data.print
b.treelist <- lapply(colnames(labiodental.data.print)[4:13], function(t) { reverse.time(b.treedata.plot) +
geom_tiplab(aes(label=label), offset=100, size=1.8) +
    geom_tippoint(aes_string(fill=t), size=2, shape=22) +
        scale_fill_manual(values=color.scheme) +
        xlim(-8200,1400) +
        annotate('text', x=-8190,y=26,label=paste("*",t)) + 
            theme_phylo(time='grid')
            })
do.call("grid.arrange", c(b.treelist, ncol=5))

# */
#' 
#' \clearpage
#' 
#' Modeling transition rates
#' ==========================
#' 
#' We model language change as CTMCs with the two states:
#' 
#' - Y: yes, the reflex of a given correspondeme is labiodental
#' - N: no, the reflex of a given correspondeme is not labiodental
#' 
#' We then estimate the rates of change *q* between two states in a Bayesian framework [@Pagel1999The-maximum-likelihood; @Huelsenbecketal2003Stochastic; @Pageletal2004Bayesian; @Ronquist2004Bayesian]:
#' 
#'  - `qYN`: instantaneous rate of loosing a labiodental realization (Yes>No)
#'  - `qNY`: instantaneous rate of gaining a labiodental realization (Yes>No)
#' 
#' In order to assess whether there is a diachronic bias in either direction, e.g. whether gains are favored over losses, we compare two models for each correspondeme:
#'
#'  - `ARD`: all rates different, i.e. rates of change towards vs away from a labiodental realization are different: `qYN`$\neq$`qNY`. This model therefore has two parameters.
#'  - `ER`: equal rates, i.e. rates of change towards vs away from a labiodental realization are the same: `qYN`=`qNY`. This model therefore has one parameter.
#' 
#' In each case, the models are compared with Bayes Factors (BF), which are defined as double the difference between the log marginal likelihood of the more complex model, i.e. the ARD model, and the log marginal likelihood of the simpler model, i.e. the ER model. BF > 2 are considered to give "positive evidence" for the better-fitting model; between 5 and 10 "strong evidence" and > 10 "very strong evidence".
#' 
#' We replicate all MCMC analyses several times[@Pageletal2014BayesTraits]  and report the median and the median absolute deviation across these replications. The parameters of each run of the analysis were determined  on the basis of visual inspections of trace plots of the likelihoods and (in the case of posterior tree sample) of the tree mixing.
#' 
#' 
#' Procedures
#' ---------
#' 
#' 
#' We use wrapper scripts which run `BayesTraitsV2` [@Pageletal2014BayesTraits] in the background and allow various steps of post-processing the output:

fit.mcmc <- function(variable, tree=ie.sum.tree, data=labiodental.data,
				prior="uniform 0 .005", assumed.mrca=NULL,
				iterations=mcmc.iterations, burn.in=mcmc.burn.in, 
				sampling.rate=mcmc.sampling.rate, covarion=F, reuse.tree=T) {	
	test.ard <- Discrete(tree, data[,c('language', variable)], "Bayesian",
		 	pa=prior, it=iterations, bi=burn.in, sa=sampling.rate,
			rm=T, fo=assumed.mrca, cv=covarion, reuse.tree=reuse.tree)
	test.er <-  Discrete(tree, data[,c('language', variable)], "Bayesian",
		 	resall='qNY', it=iterations, pa=prior, bi=burn.in, sa=sampling.rate,
			fo=assumed.mrca, cv=covarion, reuse.tree=reuse.tree)

	summary.er <- tail(test.er[,5],1)
	summary.ard <- tail(test.ard[,grep('q', colnames(test.ard), value=T)],1)
	bfres <- bftest(test.er, test.ard)
	bfres$Better <- ifelse(bfres$BetterModel %in% 'Model 1', 'ER', 'ARD')

	return(list(ARD=test.ard,
		    ER=test.er,
		    ER.vs.ARD=bfres[,c(1,3)],
		    ER.rate=summary.er,
		    ARD.rates=summary.ard
		    ))
}

fit.ml <- function(variable, tree=ie.sum.tree, data=labiodental.data, assumed.mrca=NULL) {
	# make sure we don't use the wrong trees in case they are recycled in the mcmc runs:
	if (file.exists("BT.current.tree.nex")) { file.remove('BT.current.tree.nex') } 
	test.ard <- Discrete(tree, data[,c('language', variable)], "ML",
		 	mlt=100, fo=assumed.mrca, reuse.tree=F)
	test.er <-  Discrete(tree, data[,c('language', variable)], "ML",
			 resall='qNY', mlt=100, fo=assumed.mrca, reuse.tree=F)
	return(list(ER=test.er, ARD=test.ard))
}

#' Function for tracing the MCMC iterations:
trace.plots <- function(mcmc.res, model=c('ARD','ER'), which=c('Likelihood', 'Rates', 'Trees')) {
	switch(model,
		ARD = {df <- do.call(rbind, lapply(seq_along(mcmc.res), function(r) {
			   x <- mcmc.res[[r]][['ARD']]
			   x$replicate <- r
			   return(x)
			   }))
			  },
		ER  = {df <- do.call(rbind, lapply(seq_along(mcmc.res), function(r) {
			   x <- mcmc.res[[r]][['ER']]
			   x$replicate <- r
			   return(x)
			   }))
			  }	
	)
	switch(which, 	
	Likelihood = { ggplot(df, aes(x=Iteration,y=Harmonic.Mean)) + geom_line(size=.2) + 
					facet_wrap(~replicate) + theme_bw(base_size = 12) },
	Rates =      { gather(df, rate, value, qNY:qYN) %>%
		                    ggplot(aes(x=Iteration, y=value, colour=rate)) +
		      		    	geom_line(size=.2) + theme(legend.position='top') +
					 	    guides(colour = guide_legend(override.aes = list(size=2))) +
							facet_wrap(~replicate) + theme_bw(base_size = 12)
						 },
	Trees = 	{  #ggplot(df, aes(x=Iteration, y=Tree.No)) + geom_line(size=.2)
				   ggplot(df, aes(x=Tree.No)) + geom_density(size=.2) +
			   	facet_wrap(~replicate) + theme_bw(base_size = 12) }					
	)
}


#' Some helpers for exploring and plotting the results:

model.fits <- function(m) { do.call(rbind,
	lapply(m, function(x) { data.frame("["(x,"ER.vs.ARD")) } )
	)}

fossil.fits <- function(m1, m2, model=c('ER','ARD')) { do.call(rbind,
        mapply(FUN=function(x, y) { bftest(x[[model]], y[[model]]) }, m1, m2, SIMPLIFY=F) 
        )}

rates <- function(m, model=c('ARD','ER')) {
	if(model %in% 'ARD') {mod <- paste0(model,".rates")} else {mod <- paste0(model, '.rate')}
	do.call(rbind, lapply(m, function(x) { data.frame("["(x, mod)) } )
	)}

ml.rates <- function(m) { 
	x <- do.call(rbind, sapply(m, rbind))
	x$model <- rep(c('ER', 'ARD'), nrow(x)/2)
	return(x[,c('model',grep('q', colnames(x), value=T))])
	}

transition.matrix <- function(m, aggr=c("median", "mad"), model=NULL, prob=T, years=100) {
	if(is.null(model)) {res <- sapply(m, function(x) unlist(tail(x,1)))} else {
		test_that('If you request a specific model, there should be a choice', {
			expect_true(any(names(m[[1]]) %in% c('ARD', 'ER')))
				})
	switch(model,
		ARD = {res <- sapply(m, function(x) unlist(tail(x$ARD,1)))},
		ER  = {res <- sapply(m, function(x) unlist(tail(x$ER,1)))}
		)
	}	
    switch(aggr,
		median = {rates <- apply(res[grepl('q', rownames(res)),], 1, median)},
		mad    = {rates <- apply(res[grepl('q', rownames(res)),], 1, mad)}
		) 
    # extract names, and do it in the right order so we can use it below for labeling:
	symbols <- unique(substr(gsub('q','',names(rates)),1,1))
    q.mat.diagonal <- sapply(split(rates, ceiling(seq_along(rates)/(length(symbols)-1))), 						function(x) -sum(x))  
    q.mat <- matrix(NA, ncol=length(symbols), nrow=length(symbols))
    diag(q.mat) <- q.mat.diagonal
    idx <- which(is.na(q.mat), arr.ind=T)
    q.mat[idx[order(idx[,1]),]] <- rates # reordered indices by rows first
	# perform matrix exponentiation if requested:	
	if(prob) { q.mat <- matexpo(years*q.mat) }
    dimnames(q.mat) <- list(symbols,symbols) 
	return(q.mat)
	}

plot.matrix <- function(m, label='x', legend=T, colors=color.scheme[-1],
                        node.size=20, padding=.8, weight.range=c(.5,1),
                        label.dodge=-.1, label.size=3) {
	gg <- graph_from_adjacency_matrix(adjmatrix=m, weighted=TRUE, mode="directed", diag=T)
	myarrow <- arrow(angle = 30, length = unit(0.3, "cm"), ends = "last", type = "closed")
    if(ncol(m)==2) { directions=rep(c(180,0),2) } else
			            { directions=rep(c(30,270,150),3) } # unelegant... but how else?
	p <- ggraph(gg, layout = 'linear', circular=ncol(m)>2) + 
		 geom_edge_fan(aes(width=round(weight,2), label=round(weight,2)), n=200,
			 		angle_calc='along', label_dodge=unit(label.dodge, "cm"),
                    label_size=label.size,
		 			arrow=myarrow , end_cap = circle(padding, 'cm'), 
                    start_cap = circle(padding, 'cm')
				      ) + 	
		 geom_edge_loop(aes(width=round(weight,2), direction=directions, label=round(weight,2)),arrow=myarrow,label_dodge=unit(label.dodge, "cm"), angle_calc='along',
                      label_size=label.size,
	 					end_cap = circle(padding, 'cm'), 
                        start_cap = circle(padding, 'cm')
	 		              ) +
         geom_node_point(aes(color=name), size=node.size) +                          
		 theme_graph(base_family=my.font) + theme(plot.title=element_text(hjust=.5),
         plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
	if(legend) { p + scale_edge_width(name='transitions:', range=weight.range) +
					 scale_color_manual(name=paste0("Reflex of *", label, "\nis labiodental:"), 
										values=colors) +
					 guides(color = guide_legend(override.aes = list(size = 3)))
				} else {
				 p + scale_edge_width(range=weight.range, guide=F) +
					 scale_color_manual(values=colors, guide=F) +
					 ggtitle(label)
				}			
	}

#' Functions for estimating ARD models and then using the resulting transition rates for Stochastic Character Mappping (see the explanation in Section 3.3 below):

years.per.change <- function(variable, trees=ie.trees, sum.tree=ie.sum.tree, data=labiodental.data,
                             tree.simulations=nsim, rate.replications=replications,
                             prior='uniform 0 .0025', reuse.tree=T) { 

  # fit model with BayesTraits:
  model <- replicate(rate.replications, 
	  Discrete(trees, data[,c('language', variable)],
               mode="Bayesian", pa=prior, it=mcmc.iterations, bi=mcmc.burn.in, 
               sa=mcmc.sampling.rate, cv=F, reuse.tree=reuse.tree), 
			   simplify=F)
  
  # prepare data for simmap:
  nas <- c('-','?') # catch NAs
  states <- data[, variable]
  names(states) <- data$language
  state.symbols <- setdiff(as.character(unique(states)),nas)
  n.states <- length(state.symbols)
 
  if(any(states %in% nas | is.na(states))) { # see ?make.simmap
    states <- to.matrix(states, state.symbols)
	states[rowSums(states)==0,] <- rep(1/n.states) # equal probability
  } 
  
  # extract Q matrix:  
  q.mat <- transition.matrix(model, 'median', prob=F)
  
  # estimate simmaps with a flat root prior:
  treemaps <- make.simmap(sum.tree, states, Q=q.mat, pi='estimated', 
                          nsim=tree.simulations, message=F)

  # extract number of changes, translate into how many years a change takes:
  simmap.description <- describe.simmap(treemaps, message=F)
  median.n.changes <- median(simmap.description$count[,'N']) # N = N(transitions)
  min.n.changes <- min(simmap.description$count[,'N'])
  max.n.changes <- max(simmap.description$count[,'N'])
  return(list(simmap.description=simmap.description,
	  		  median.years = round(sum(sum.tree$edge.length)/median.n.changes),
              min.years = round(sum(sum.tree$edge.length)/max.n.changes),
			  max.years = round(sum(sum.tree$edge.length)/min.n.changes)
			  ))
}
	
est.changes <- function(simmap.description, sum.tree) {
	apply(sapply(simmap.description, function(x) x$count[,'N']), 2, function(col) { 
        sum(sum.tree$edge.length)/col
		}) %>% 	as.data.frame %>% 
	gather(correspondeme, years)
	}
		
	
#' 
#' Parameters
#' ----------
#' 
#' *Replications:*
B=10

#' *Replications of covarion models*:
B.cv=3

#' *MCMC iterations:*
mcmc.iterations=1e+07

#' *Burn in:*
mcmc.burn.in=1e+05

#' *Sampling rate:*
mcmc.sampling.rate=1e+03


#+ stepcheck1, include=F
run=T # set to F so that the next steps won't be executed

#' 
#' 
#' Evaluation of priors
#' --------------------
#' 
#' We defined the prior distribution for each model so that the estimated number of changes is within a realistic ballpark. To assess this ballpark we explored transition rate estimates obtained with priors of $\mathcal{U}(0,.0025)$, $\mathcal{U}(0,.005)$ and $\mathcal{U}(0,.01)$ by performing *Stochastic Character Mapping* (see Section 5 below). For this purpose we focus on models allowing for rate differences (ARD). This allows a larger range of rates to explore. We then also compare the rates obtained with these priors with rates estimated from Maximum-Likelihood models.
#' 
#' **Parameters:**
#' 
#' *Stochastic character mappings:*
nsim=100

#' *Replications of rate estimates:*
replications=5

#' 
#' ### Speed of change allowed under various prior distributions on rates
#' 
correspondemes <- c('bh','p','w','m','gw','gwh','kw','dh','ky','b')
#' 
#' 
#' 
#' **Phylogenies from Chang et al.**
#' 
#+ priors-c, cache=T, results='hide', echo=T, eval=run
ypc_0_.0025c <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.c.trees, sum.tree=ie.c.sum.tree, 
						prior='uniform 0 .0025') })
ypc_0_.005c <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.c.trees, sum.tree=ie.c.sum.tree, 
						prior='uniform 0 .005') })
ypc_0_.01c <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.c.trees, sum.tree=ie.c.sum.tree, 
						prior='uniform 0 .01') })
#' 
#' 
#' **Phylogenies from Bouckaert et al.**
#' 
#+ priors-b, cache=T, results='hide', echo=T, eval=run
ypc_0_.0025b <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.b.trees, sum.tree=ie.b.sum.tree, 
						prior='uniform 0 .0025') })
ypc_0_.005b <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.b.trees, sum.tree=ie.b.sum.tree, 
						prior='uniform 0 .005') })
ypc_0_.01b <- sapply(correspondemes, function(crsp) {
					years.per.change(crsp, trees=ie.b.trees, sum.tree=ie.b.sum.tree, 
						prior='uniform 0 .01') })

#+ priors-b-sum, fig.cap='Estimated number of years per change per prior assumptions about rates of change', eval=run, echo=F, fig.pos='b!'
pr.c <- gdata::combine(est.changes(ypc_0_.0025b[1,], ie.b.sum.tree), 
			   est.changes(ypc_0_.005b[1,],  ie.b.sum.tree), 
			   est.changes(ypc_0_.01b[1,],   ie.b.sum.tree)) %>% 
	mutate(prior=gsub('(.*)(\\.\\d{1,4})([bc].*)', 'U(0, \\2)', source),
           correspondeme=factor(correspondeme, levels=correspondemes)) 
pr.b <- gdata::combine(est.changes(ypc_0_.0025c[1,], ie.c.sum.tree), 
			   est.changes(ypc_0_.005c[1,],  ie.c.sum.tree), 
			   est.changes(ypc_0_.01c[1,],   ie.c.sum.tree)) %>% 
	mutate(prior=gsub('(.*)(\\.\\d{1,4})([bc].*)', 'U(0, \\2)', source),
           correspondeme=factor(correspondeme, levels=correspondemes))           

gdata::combine(pr.b,pr.c) %>%       
ggplot(aes(x=correspondeme, y=years, fill=source.1)) + 
	geom_boxplot() +
    # annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=30, fill='red', alpha=.2) +
	facet_wrap(~prior) +
	scale_y_log10(breaks=c(10, 100, 1000, 10000, 100000),
         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(name='',
                    labels=c("bʰ","p","w","m","gʷ","gʷʰ","kʷ","dʰ","kʲ","b")) +
    scale_fill_manual(name='phylogeny:', values=c('grey','grey95'), 
                      labels=c('Bouckaert et al', 'Chang et al.')) +
    theme(legend.position='bottom')                                 
    
#' 
#' ### Maximum Likelihood estimates of rates
#'   
#' 
#' **Phylogenies from Chang et al.**
#'   
#+ ml-c, cache=T, eval=run
(ml.estimates.c <- sapply(correspondemes, function(crsp) {
		ml <- replicate(B, fit.ml(crsp, tree=ie.c.sum.tree), simplify=F)
		rates <- unlist(ml.rates(ml)[,-1])
		return(data.frame(median=round(median(rates),5), max=round(max(rates),5)))
	 	}))
#' 
#' **Phylogenies from Bouckaert et al.**
#' 
#+ ml-b, cache=T, eval=run
(ml.estimates.b <- sapply(correspondemes, function(crsp) {
		ml <- replicate(B, fit.ml(crsp, tree=ie.b.sum.tree), simplify=F)
		rates <- unlist(ml.rates(ml)[,-1])
		return(data.frame(median=round(median(rates),5), max=round(max(rates),5)))
	 	}))


#+ stepcheck1b, include=F
run=T # set to F so that the next steps won't be executed

#' 		
#' Individual analyses
#' -------------------
#' 
#' 
#' Given the results of the previous section, a prior of $\mathcal{U}(0,.01)$ seems reasonable, i.e. it is large enough without allowing unrealistically fast or unrealistically slow rates of change:
prior="uniform 0 .01"
#' 
#' 
#' ### Phylogenies from Chang et al.
#' 
#' Reflexes of *$b^h$:
#' 
#+ mcmc-bh-c, eval=run, cache=T
bh.mcmc.c <- replicate(B, fit.mcmc('bh',   tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-bh-res-c, eval=run, include=F
# trace.plots(bh.mcmc.c, model='ER', which='Likelihood')
# trace.plots(bh.mcmc.c, model='ARD', which='Likelihood')
model.fits(bh.mcmc.c)

#' 
#'  Reflexes of *$p$:
#' 
#+ mcmc-p-c, eval=run, cache=T
p.mcmc.c <- replicate(B, fit.mcmc('p',     tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-p-res-c, eval=run, include=F
# trace.plots(p.mcmc.c, model='ER', which='Likelihood')
# trace.plots(p.mcmc.c, model='ARD', which='Likelihood')
model.fits(p.mcmc.c) 

# /*
# plot.matrix(transition.matrix(p.mcmc.c, 'median', 'ER', years=1000), label='p')
# */

#' 
#'  Reflexes of *$w$:
#' 
#+ mcmc-w-c, eval=run, cache=T
w.mcmc.c <- replicate(B, fit.mcmc('w',     tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-w-res-c, eval=run, include=F
# trace.plots(w.mcmc.c, model='ER', which='Likelihood')
# trace.plots(w.mcmc.c, model='ARD', which='Likelihood')
model.fits(w.mcmc.c) 

# /*
# plot.matrix(transition.matrix(w.mcmc.c, 'median', 'ER', years=1000), label='w')
# */

#' 
#'  Reflexes of *$m$:
#' 
#+ mcmc-m-c, eval=run, cache=T
m.mcmc.c <- replicate(B, fit.mcmc('m',     tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-m-res-c, eval=run, include=F
# trace.plots(m.mcmc.c, model='ER', which='Likelihood')
# trace.plots(m.mcmc.c, model='ARD', which='Likelihood')
model.fits(m.mcmc.c) 

#' 
#'  Reflexes of *$g^w$:
#' 
#+ mcmc-g-c, eval=run, cache=T
gw.mcmc.c <- replicate(B, fit.mcmc('gw',   tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-gw-res-c, eval=run, include=F
# trace.plots(gw.mcmc.c, model='ER', which='Likelihood')
# trace.plots(gw.mcmc.c, model='ARD', which='Likelihood')
model.fits(gw.mcmc.c) 

# /*
# plot.matrix(transition.matrix(gw.mcmc, 'median', 'ER', years=1000), label='gʷ')
# */

#' 
#'  Reflexes of *$g^{wh}$:
#' 
#+ mcmc-gwh-c, eval=run, cache=T
gwh.mcmc.c <- replicate(B, fit.mcmc('gwh', tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-gwh-res-c, eval=run, include=F
# trace.plots(gwh.mcmc.c, model='ER', which='Likelihood')
# trace.plots(gwh.mcmc.c, model='ARD', which='Likelihood')
model.fits(gwh.mcmc.c)

# /*
# plot.matrix(transition.matrix(gwh.mcmc, 'median', 'ER', years=1000), label='gʷʰ')
# */

#' 
#'  Reflexes of *$k^w$:
#' 
#+ mcmc-kw-c, eval=run, cache=T
kw.mcmc.c <- replicate(B, fit.mcmc('kw',   tree=ie.c.trees, prior=prior), simplify=F)

#+ mcmc-kw-res-c, eval=run, include=F
# trace.plots(kw.mcmc.c, model='ER', which='Likelihood')
# trace.plots(kw.mcmc.c, model='ARD', which='Likelihood')
model.fits(kw.mcmc.c)

#' 
#'  Reflexes of *$d^h$:
#' 
#+ mcmc-dh-c, eval=run, cache=T
dh.mcmc.c <- replicate(B, fit.mcmc('dh',   tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-dh-res-c, eval=run, include=F
# trace.plots(dh.mcmc.c, model='ER', which='Likelihood')
# trace.plots(dh.mcmc.c, model='ARD', which='Likelihood')
model.fits(dh.mcmc.c)
#' 
#'  Reflexes of *$k^j$:
#' 
#+ mcmc-k-c, eval=run, cache=T
ky.mcmc.c <- replicate(B, fit.mcmc('ky',     tree=ie.c.trees, prior=prior), simplify=F)

#+ mcmc-k-res-c, eval=run, include=F
# trace.plots(k.mcmc.c, model='ER', which='Likelihood')
# trace.plots(k.mcmc.c, model='ARD', which='Likelihood')
model.fits(ky.mcmc.c)

#' 
#'  Reflexes of *$b$:
#' 
#+ mcmc-b-c, eval=run, cache=T
b.mcmc.c <- replicate(B, fit.mcmc('b',     tree=ie.c.trees, prior=prior), simplify=F)
#+ mcmc-b-res-c, eval=run, include=F
# trace.plots(b.mcmc.c, model='ER', which='Likelihood')
# trace.plots(b.mcmc.c, model='ARD', which='Likelihood')
model.fits(b.mcmc.c) 

#' 
#' 
#' ### Phylogenies from Bouckaert et al.
#' 
#+ echo=F, results='hide'
# Make sure we use the right trees:
if (file.exists("BT.current.tree.nex")) { file.remove('BT.current.tree.nex') } 
#' 
#'  Reflexes of *$b^h$:
#' 
#+ mcmc-bh-b, eval=run, cache=T
bh.mcmc.b <- replicate(B, fit.mcmc('bh',   tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-bh-res-b, eval=run, include=F
# trace.plots(bh.mcmc.b, model='ER', which='Likelihood')
# trace.plots(bh.mcmc.b, model='ARD', which='Likelihood')
model.fits(bh.mcmc.b)

#' 
#'  Reflexes of *$p$:
#' 
#+ mcmc-p-b, eval=run, cache=T
p.mcmc.b <- replicate(B, fit.mcmc('p',     tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-p-res-b, eval=run, include=F
# trace.plots(p.mcmc.b, model='ER', which='Likelihood')
# trace.plots(p.mcmc.b, model='ARD', which='Likelihood')
model.fits(p.mcmc.b) 

#' 
#'  Reflexes of *$w$:
#' 
#+ mcmc-w-b, eval=run, cache=T
w.mcmc.b <- replicate(B, fit.mcmc('w',     tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-w-res-b, eval=run, include=F
# trace.plots(w.mcmc.b, model='ER', which='Likelihood')
# trace.plots(w.mcmc.b, model='ARD', which='Likelihood')
model.fits(w.mcmc.b) 

#' 
#'  Reflexes of *$m$:
#' 
#+ mcmc-m-b, eval=run, cache=T
m.mcmc.b <- replicate(B, fit.mcmc('m',     tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-m-res-b, eval=run, include=F
# trace.plots(m.mcmc.b, model='ER', which='Likelihood')
# trace.plots(m.mcmc.b, model='ARD', which='Likelihood')
model.fits(m.mcmc.b) 

#' 
#'  Reflexes of *$g^w$:
#' 
#+ mcmc-g-b, eval=run, cache=T
gw.mcmc.b <- replicate(B, fit.mcmc('gw',   tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-gw-res-b, eval=run, include=F
# trace.plots(gw.mcmc.b, model='ER', which='Likelihood')
# trace.plots(gw.mcmc.b, model='ARD', which='Likelihood')
model.fits(gw.mcmc.b) 

#' 
#'  Reflexes of *$g^{wh}$:
#' 
#+ mcmc-gwh-b, eval=run, cache=T
gwh.mcmc.b <- replicate(B, fit.mcmc('gwh', tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-gwh-res-b, eval=run, include=F
# trace.plots(gwh.mcmc.b, model='ER', which='Likelihood')
# trace.plots(gwh.mcmc.b, model='ARD', which='Likelihood')
model.fits(gwh.mcmc.b)

#' 
#'  Reflexes of *$k^w$:
#' 
#+ mcmc-kw-b, eval=run, cache=T
kw.mcmc.b <- replicate(B, fit.mcmc('kw',   tree=ie.b.trees, prior=prior), simplify=F)

#+ mcmc-kw-res-b, eval=run, include=F
# trace.plots(kw.mcmc.b, model='ER', which='Likelihood')
# trace.plots(kw.mcmc.b, model='ARD', which='Likelihood')
model.fits(kw.mcmc.b)

#' 
#'  Reflexes of *$d^h$:
#' 
#+ mcmc-dh-b, eval=run, cache=T
dh.mcmc.b <- replicate(B, fit.mcmc('dh',   tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-dh-res-b, eval=run, include=F
# trace.plots(dh.mcmc.b, model='ER', which='Likelihood')
# trace.plots(dh.mcmc.b, model='ARD', which='Likelihood')
model.fits(dh.mcmc.b)

#' 
#'  Reflexes of *$k^j$:
#' 
#+ mcmc-k-b, eval=run, cache=T
ky.mcmc.b <- replicate(B, fit.mcmc('ky',     tree=ie.b.trees, prior=prior), simplify=F)

#+ mcmc-k-res-b, eval=run, include=F
# trace.plots(k.mcmc.b, model='ER', which='Likelihood')
# trace.plots(k.mcmc.b, model='ARD', which='Likelihood')
model.fits(ky.mcmc.b)

#' 
#'  Reflexes of *$b$:
#' 
#+ mcmc-b-b, eval=run, cache=T
b.mcmc.b <- replicate(B, fit.mcmc('b',     tree=ie.b.trees, prior=prior), simplify=F)
#+ mcmc-b-res-b, eval=run, include=F
# trace.plots(b.mcmc.b, model='ER', which='Likelihood')
# trace.plots(b.mcmc.b, model='ARD', which='Likelihood')
model.fits(b.mcmc.b) 

#' 
#' Results
#' -------
#' 
#' Figure S4.4 summarizes the models in terms of the logarithm of the Bayes Factor (logBF) across replications and phylogenies. logBF is defined as double the difference between the log marginal likelihood of the more complex model (assuming unequal rates) and the log marginal likelihood of the simpler model (assuming equal rates). Values above 2 are usually considered to support the complex model and this is indicated by light blue shading in the graphs. Only the correspondemes $\text{*}m$, $\text{*}g^w$, and $\text{*}b$ show evidence for more complex models with unequal rates. The rate estimates themselves are plotted in Figure S4.5. The three correspondemes with unequal rates are more likely to lose labiodental reflexes than to gain them.
#' 
#' 
#+ chains, echo=F, include=F, eval=F

names(correspondemes) <- c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")
chain.list.c <- lapply(seq_along(correspondemes), function(i) {
    correspondeme <- correspondemes[i]
    print(correspondeme)
    model <- get(paste(tolower(correspondeme), 'mcmc.c', sep="."))
    tab <- table(model.fits(model)$ER.vs.ARD.Better)
    best <- names(tab[tab==max(tab)])
    print(best)
    print(2*(tail(model$ARD$Harmonic.Mean,1)-tail(model$ER$Harmonic.Mean,1)))
    trans.m <- transition.matrix(model, 'median', best, years=20000)
    print(round(trans.m,2))
    m <- plot.matrix(trans.m, legend=F, label=paste(names(correspondemes[i])),
        node.size=10, padding=.4, weight.range=c(.5,1), label.dodge=0.2)
})
quartz(height=9, width=15)
do.call("grid.arrange", c(chain.list.c, ncol=5))

#' 
#+ sum-ctmc, echo=F, eval=run, include=F
# kable(data.frame(Correspondeme=c('$\\text{*}b^h$','$\\text{*}p$','$\\text{*}w$','$\\text{*}m$','$\\text{*}g^w$','$\\text{*}g^{wh}$','$\\text{*}k^w$','$\\text{*}d^h$','$\\text{*}k^j$','$\\text{*}b$'),
#    C.best=sapply(correspondemes, function(t) {
#        paste(sort(unique(model.fits(get(paste(tolower(t), 'mcmc.c', sep=".")))[,2])),
#                collapse=', ') }),
#    C.median=sapply(correspondemes,
#            function(t) median(model.fits(get(paste(tolower(t), 'mcmc.c', sep=".")))[,1])),
#    C.mad=sapply(correspondemes,
#            function(t) mad(model.fits(get(paste(tolower(t), 'mcmc.c', sep=".")))[,1])),
#    B.best=sapply(correspondemes, function(t) {
#        paste(sort(unique(model.fits(get(paste(tolower(t), 'mcmc.b', sep=".")))[,2])),
#                collapse=', ') }),
#    B.median=sapply(correspondemes,
#            function(t) median(model.fits(get(paste(tolower(t), 'mcmc.b', sep=".")))[,1])),
#    B.mad=sapply(correspondemes,
#            function(t) mad(model.fits(get(paste(tolower(t), 'mcmc.b', sep=".")))[,1]))
#            ), digits=2, row.names=F,
#            caption="Summary of models across correspondemes and trees")
#' 
#'   
#+ sum-ctmc-plot-c, echo=F, eval=run, fig.cap='Bayes Factor evidence for different-rate model compared to equal-rate model. The blue shading covers the range conventionally considered not to give evidence.' 
bf.c <- sapply(correspondemes, function(t) {
    sapply(get(paste(tolower(t), 'mcmc.c', sep=".")), function(m) {            
        2*(tail(m$ARD$Harmonic.Mean,1)-tail(m$ER$Harmonic.Mean,1))
        })
    }) %>% as.data.frame %>% 
    gather(type, logBF)  %>% 
    mutate(correspondeme=factor(type, levels=correspondemes,
            labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b"))) 
bf.b <- sapply(correspondemes, function(t) {
    sapply(get(paste(tolower(t), 'mcmc.b', sep=".")), function(m) {            
        2*(tail(m$ARD$Harmonic.Mean,1)-tail(m$ER$Harmonic.Mean,1))
        })
    }) %>% as.data.frame %>% 
    gather(type, logBF)  %>% 
    mutate(correspondeme=factor(type, levels=correspondemes,
            labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")))            

gdata::combine(bf.b,bf.c) %>% 
ggplot(aes(x=correspondeme, y=logBF, fill=source)) + geom_boxplot() +
    scale_y_continuous(breaks=c(-10,-5,-2,0,2,5,10)) + 
    scale_fill_manual(name='phylogeny:', values=c('grey','grey95'), 
                                         labels=c('Bouckaert et al','Chang et al.')) + 
    xlab('') +
    annotate('rect',xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=2, fill='lightblue', alpha=.20)

#' 
#+ rates-ctmc-plot, echo=F, eval=run, fig.cap='Estimated rates of gain and loss of labiodentals in the best models. Gray colors indicate that the best model assumes equal rates of gain and loss. Dark colors indicate phylogenies from Bouckaert et al (B), light colors phylogenies from Chang et al. (C).' 

results.df <- data.frame(best=sapply(correspondemes, function(t) { 
    	   paste(sort(unique(model.fits(get(paste(tolower(t), 'mcmc.c', sep=".")))[,2])), 
    	   		collapse=', ') }))

rates.df <- mapply(function(c,m) {
                    best.model <- ifelse(m %in% 'ARD','ARD','ER')
                    sapply(get(paste(c, 'mcmc.c', sep=".")), function(x) {
                            data.frame(qNY=tail(x[[best.model]]$qNY,1), 
                                       qYN=tail(x[[best.model]]$qYN,1))
                        })
                    }, 
               rownames(results.df), results.df$best, SIMPLIFY=F) %>% 
    as.data.frame %>% 
    gather(type, value)  %>% 
    separate(type, into=c('correspondeme', 'run'), sep='\\.') %>% 
    mutate(value=unlist(value),
            rate=rep(c('gain','loss'), 10*10), 
            correspondeme=factor(correspondeme, levels=correspondemes,
               labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")))

results.b.df <- data.frame(best=sapply(correspondemes, function(t) { 
    	   paste(sort(unique(model.fits(get(paste(tolower(t), 'mcmc.b', sep=".")))[,2])), 
    	   		collapse=', ') }))

rates.b.df <- mapply(function(c,m) {
                    best.model <- ifelse(m %in% 'ARD','ARD','ER')
                    sapply(get(paste(c, 'mcmc.b', sep=".")), function(x) {
                            data.frame(qNY=tail(x[[best.model]]$qNY,1), 
                                       qYN=tail(x[[best.model]]$qYN,1))
                        })
                    }, 
               rownames(results.b.df), results.b.df$best, SIMPLIFY=F) %>% 
    as.data.frame %>% 
    gather(type, value)  %>% 
    separate(type, into=c('correspondeme', 'run'), sep='\\.') %>% 
    mutate(value=unlist(value),
            rate=rep(c('gain','loss'), 10*10), 
            correspondeme=factor(correspondeme, levels=correspondemes,
               labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")))

rates.all.df <- gdata::combine(rates.df, rates.b.df) %>% mutate(rate.source=paste(rate,source,'.'))

ggplot(rates.all.df, aes(x=correspondeme, y=value)) + 
    geom_point(size=0, alpha=0) + # just a dummy to ensure consistent level ordering
    geom_boxplot(data=subset(rates.all.df, correspondeme %in% c('*m','*b', '*gʷ')), aes(fill=rate.source)) +
    geom_boxplot(data=subset(rates.all.df, !correspondeme %in% c('*m','*b', '*gʷ')), aes(fill=source)) +
    ylab('rate estimate') + xlab('') +
    scale_fill_manual(name='', values=c('red4','red1', 'blue4', 'blue1', 'grey', 'grey95'),
              labels=c('gain (B)', 'gain (C)', 'loss (B)', 'loss (C)', 'rate (B)', 'rate (C)')
                               )


#+ gw, echo=F, eval=run
gw.gain <- subset(rates.b.df, correspondeme %in% '*gʷ' & rate %in% 'gain')
gw.loss <- subset(rates.b.df, correspondeme %in% '*gʷ' & rate %in% 'loss')

#' 
#' 
#' \clearpage
#' 
#' Ancestral state estimates
#' ========================
#' 
#' CTMC models include an estimate on the probability distribution of the root value. We extract this from the best-fitting models and report the results in Figure S4.6.
#' 
#+ ace-df, cache=T
posterior.size=(mcmc.iterations-mcmc.burn.in)/mcmc.sampling.rate
ace.df <- data.frame(
           phylogeny=c(rep("Chang et al.'s phylogeny", 
                             B*length(correspondemes)*posterior.size), 
                       rep("Bouckaert et al.'s phylogeny", 
                             B*length(correspondemes)*posterior.size)),
           correspondeme=factor(rep(
                               unlist(lapply(correspondemes,rep,B*posterior.size)),
                               2), 
                               levels=correspondemes,
                               labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")),
           Prob=c(unlist(lapply(correspondemes, function(crsp) {
                   sapply(get(paste0(crsp, '.mcmc.c')), function(x) {
                               if(crsp %in% c('*m','*b','*gʷ')) { x$ARD$Root.P.Y. } else {
                               x$ER$Root.P.Y.
                               }
                           }
                        ) })),
                   unlist(lapply(correspondemes, function(crsp) {
                    sapply(get(paste0(crsp, '.mcmc.b')), function(x) {
                               if(crsp %in% c('*m','*b','*gʷ')) { x$ARD$Root.P.Y. } else {
                               x$ER$Root.P.Y.
                               }
                           }
                        ) }))
                      )
            ) 

#' Figure S4.6 suggests that for 8 correspondemes, the probability of a labiodental realization in proto-Indo-European was very low. For two correspondemes --- $\text{*}w$ and $\text{*}b$ --- the ancestral state cannot be sufficiently resolved and is estimated at or close to 50\% probabilities. In these two cases, there is no decisive evidence in favor of nor against a labiodental realization in proto-Indo-European. One way of resolving this issue is to fit and compare models with vs without assumed proto-states.
#' 
#+ ace-plot, fig.cap='Posterior probabilities of labiodentals in PIE, based on the best-fitting models', echo=F, fig.pos='!b'          
ggplot(ace.df, aes(x=Prob, y=correspondeme)) + facet_wrap(~phylogeny) +
    geom_density_ridges(scale=.99, rel_min_height = 0.01, alpha=.5, size=.3) +
    xlab('') +
    ylab('') 

# /*
# posterior.size=(mcmc.iterations-mcmc.burn.in)/mcmc.sampling.rate
# ace.df <- data.frame(
#            correspondeme=factor(unlist(lapply(correspondemes,rep,B*posterior.size)),
#                                levels=correspondemes,
#                                labels=c("*bʰ","*p","*w","*m","*gʷ","*gʷʰ","*kʷ","*dʰ","*kʲ","*b")),
#            Prob=unlist(lapply(correspondemes, function(crsp) {
#                    sapply(get(paste0(crsp, '.mcmc.c')), function(x) {
#                                if(crsp %in% c('*m','*b','*gʷ')) { x$ARD$Root.P.Y. } else {
#                                x$ER$Root.P.Y.
#                                }
#                            }
#                         ) }))
#             )
#
# ggplot(arrange(ace.df, desc(correspondeme)), aes(x=Prob, y=correspondeme)) +
#     geom_density_ridges(
#         scale=.99, rel_min_height = 0.01,
#         fill='coral1', color='black', size=.3
#         ) +
#         xlab('Posterior probability of proto-labiodental') +
#         ylab('') +
#         theme(panel.background = element_rect(fill = 'ivory1', colour = 'white'),
#               panel.grid.major = element_line(colour = "black", size=.1),
#               axis.line=element_line(color='black'),
#               axis.text.y = element_text(color='black',size=14),
#               axis.text.x=element_text(size=14),
#               axis.title=element_text(size=18))
# */
#' 
#' \clearpage
#' 
#' Reconstruction comparing assumed proto-states
#' --------------------------------------------
#' 
#' ### Based on the phylogenies from Chang et al.
#' 
#' **Reconstruction of $\text{*}w$**:
#' 
#+ mcmc-w-f-c, eval=run, cache=T
PIE.had.labiodentals.c <- paste("Y", paste(ie.c.trees[[1]]$tip.label, collapse=" "))
PIE.lacked.labiodentals.c <- paste("N", paste(ie.c.trees[[1]]$tip.label, collapse=" "))

w.mcmc.c.Y <- replicate(B, fit.mcmc('w', tree=ie.c.trees, prior=prior, 
                            assumed.mrca=PIE.had.labiodentals.c), simplify=F)
w.mcmc.c.N <- replicate(B, fit.mcmc('w', tree=ie.c.trees, prior=prior, 
                            assumed.mrca=PIE.lacked.labiodentals.c), simplify=F)

#+ mcmc-w-res-f-c, eval=run
#' The best fitting model is:
unique(model.fits(w.mcmc.c.Y)[,2])
median(model.fits(w.mcmc.c.Y)[,1])

#' There is no evidence for or against a model that assumes labiodentals:
fossil.fits(w.mcmc.c.Y, w.mcmc.c.N, model='ER')

#' **Reconstruction of $\text{*}b$**:
#' 
#+ mcmc-b-f-c, eval=run, cache=T
b.mcmc.c.Y <- replicate(B, fit.mcmc('b', tree=ie.c.trees, prior=prior, 
                            assumed.mrca=PIE.had.labiodentals.c), simplify=F)
b.mcmc.c.N <- replicate(B, fit.mcmc('b', tree=ie.c.trees, prior=prior, 
                            assumed.mrca=PIE.lacked.labiodentals.c), simplify=F)

#+ mcmc-b-res-f-c, eval=run
#' The best fitting model is:
unique(model.fits(b.mcmc.c.Y)[,2])
median(model.fits(b.mcmc.c.Y)[,1])
#' There is no evidence for or against a model that assumes labiodentals:
fossil.fits(b.mcmc.c.Y, b.mcmc.c.N, model='ARD')

#' 
#' ### Based on the phylogenies from Bouckaert et al.
#' 
#' **Reconstruction of $\text{*}w$**:
#' 
#+ echo=F, results='hide'
# Make sure we use the right trees:
if (file.exists("BT.current.tree.nex")) { file.remove('BT.current.tree.nex') } 

#+ mcmc-w-f-b, eval=run, cache=T
PIE.had.labiodentals.b <- paste("Y", paste(ie.b.trees[[1]]$tip.label, collapse=" "))
PIE.lacked.labiodentals.b <- paste("N", paste(ie.b.trees[[1]]$tip.label, collapse=" "))

w.mcmc.b.Y <- replicate(B, fit.mcmc('w', tree=ie.b.trees, prior=prior, 
                            assumed.mrca=PIE.had.labiodentals.b), simplify=F)
w.mcmc.b.N <- replicate(B, fit.mcmc('w', tree=ie.b.trees, prior=prior, 
                            assumed.mrca=PIE.lacked.labiodentals.b), simplify=F)

#+ mcmc-w-res-f-b, eval=run
#' The best fitting model is:
unique(model.fits(w.mcmc.b.Y)[,2])
median(model.fits(w.mcmc.b.Y)[,1])
#' There is no evidence for or against a model that assumes labiodentals:
fossil.fits(w.mcmc.b.Y, w.mcmc.b.N, model='ER')

#' **Reconstruction of $\text{*}b$**:
#' 
#+ mcmc-b-f-b, eval=run, cache=T
b.mcmc.b.Y <- replicate(B, fit.mcmc('b', tree=ie.b.trees, prior=prior, 
                            assumed.mrca=PIE.had.labiodentals.b), simplify=F)
b.mcmc.b.N <- replicate(B, fit.mcmc('b', tree=ie.b.trees, prior=prior, 
                            assumed.mrca=PIE.lacked.labiodentals.b), simplify=F)

#+ mcmc-b-res-f-b, eval=run
#' The best fitting model is:
unique(model.fits(b.mcmc.b.Y)[,2])
median(model.fits(b.mcmc.b.Y)[,1])
#' There is no evidence for or against a model that assumes labiodentals:
fossil.fits(b.mcmc.b.Y, b.mcmc.b.N, model='ARD')

#'   
#' Summary
#' -------
#' 
#' In summary, there is no evidence for or against positing labiodentals as ancestral states of the $\text{*}w$ and $\text{*}b$ correspondemes. In the case of $\text{*}b$ we earlier observed a diachronic bias towards *loss* of labiodental realizations (Figure S4.5). The rate estimates for $\text{*}w$ suggest general instability, with equally high rates of gain and loss. This confirms the qualitative observation that this correspondeme varies substantially in closely related languages (e.g. German vs. English, Spanish vs. Gallo-Romance) and its realization is highly unstable in some extant languages. For example in Hindi, there is substantial variation between [w], [$\beta$], or [v] [@Shapiro2003Hindi].
#' 
#' \clearpage              
#' 
#' 
#+ stepcheck2, include=F
run=T # set to F so that the next steps won't be executed

#' Stochastic Character Mapping
#' ============================
#' 
#' Stochastic Character Mapping [@Nielsen2002Mapping; @Huelsenbecketal2003Stochastic] samples a posterior distribution of character histories given a phylogeny, the current states, and transition rates. We estimate maps using the rate estimates from the best-fitting models estimated in the previous section.  The maps are estimated using the algorithm in `phytools::make.simmap` [@Revell2012phytools] and relies on Maximum Clade Credibility (MCC) summary trees for easier visualization as probability changes along each branch. Ancestral states at nodes however are estimates using posterior phylogeny samples and are plotted as pie charts.
#' 
#' In all models, we assume constant rates of change over time. This assumption is conservative and disfavors our hypothesis because it will over-estimate the age of innovations. Non-constant rates of change would allow for a late acceleration of innovations.
#' 
#' Procedures
#' ----------

compute.simmap <- function(variable, tree=ie.c.sum.tree, data=labiodental.data,
         model=c('ER','ARD'), iterations=B, resolution=res) {
		# get the right model based on the tree
		if(grepl('c', deparse(substitute(tree)))) {type='c'} else {type='b'}	

		# get Q matrix from BayesTraits results:					 
		mcmc.results <- get(paste0(tolower(variable), '.mcmc.', type))	
		Q.matrix <- transition.matrix(mcmc.results, 'median', model, prob=F)
        # prepare data for simmap, catching NAs (which will be estimated)
		nas <- c('-','?') # catch NAs
 	   	states <- data[, variable]
 	   	names(states) <- data$language
		state.symbols <- sort(setdiff(as.character(unique(states)),nas))
		n.states <- length(state.symbols)
            # estimating NAs with a uniform prior:
  	  		if(any(states %in% nas | is.na(states))) { 
    				states <- to.matrix(states, state.symbols)
					states[rowSums(states)==0,] <- rep(1/n.states) 
  			  	} 
	
    	# make simmaps with a flat root prior:
		treemaps <- make.simmap(tree, states, model=model, Q=Q.matrix, 
                    pi='estimated', nsim=iterations, message=F)
		}

#' The following function estimates SCMs on a single tree, but it adds node estimates based on a posterior tree sample:

plot.density <- function(variable='bh', tree=ie.c.sum.tree, 
                         nodes=T, multi.ace=T, legend=2000, lwd=c(8,8)) {
        if(grepl('c', deparse(substitute(tree)))) {type='c'} else {type='b'}
        # get the density map:
        density.map <- get(paste0(tolower(variable), '.density.', type))	
		density.map$tree <- ladderize.simmap(density.map$tree) 
		plot(density.map, outline=T, lwd=lwd, 
			fsize=c(1,1), ftype='reg', legend=legend)
        if(nodes) {
            # get the nodes that are ancestral and different from ancient languages, 
            # depending on the tree:
            switch(type, 
                c = nodes <- c(53:59, 70, 79, 84, 99, 100, 103),
                 # nodes <- c(53:61, 66:73, 75, 77, 80:84, 86:95, 97:100, 103)
                b = nodes  <- c(53:59, 65, 70, 79, 84, 99, 100, 102:103) 
                )
            # reconstruct based on stochastic mapping on a single tree or on a sample of trees:
            if(multi.ace) {scm <- '.multi.simmaps.'} else {scm <- '.simmaps.'}
            simmaps <- get(paste0(tolower(variable), scm, type))	
            node.aces <- as.data.frame(describe.simmap(simmaps)$ace)
            node.aces$node <- rownames(node.aces)
            # and plot:
            nodelabels(node=nodes, pie=as.matrix(node.aces[node.aces$node %in% nodes, 1:2]),
                                   piecol=c('blue','red'), cex=.7)   
        }    
  }	

#' 
#' 
#' Parameters
#' ----------
#' 
#' *Number of simulated histories per tree:*
B=100
#' *Resolution (bins)*: the number of time intervals (bins) is fixed by the resolution parameter of `phytools::densityMap`, which we set to:
res=200
#' Setting the parameter to this value results in bins of about 30 years, roughly reflecting the state at one generation of speakers.
#' 
#' 
#' Individual analyses
#' -------------------
#' 
#' A full analysis on the posterior tree sample is computationally extremely expensive (*B* simulations $\times$ 10 correspondemes $\times$ 20,000 trees). In response to this, we use only the last 1000 trees from the posterior sample of each phylogeny: 
#+tree-samples, cache=T, eval=run
ie.c.trees.sample <- tail(ie.c.trees, 1000)
ie.b.trees.sample <- tail(ie.b.trees, 1000)

#' 
#' The results of this is displayed as pie charts at the nodes. But the stochastic character maps along the branches are based on MCC summary trees.
#' 
#' ### Reconstruction of labiodental reflexes for *$b^h$
#' 
#' With phylogenies estimated by Chang et al.:
#+ bh.simmap-c, cache=T, eval=run, results='hide'
correspondeme='bh'
bh.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
bh.density.c <- densityMap(bh.simmaps.c, plot=F, res=res)
bh.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ bh-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$b^h$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)

#+ stepcheck2b, include=F
run=T # set to F so that the next steps won't be executed

#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ bh.simmap-b, cache=T, eval=run, results='hide'
bh.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
bh.density.b <- densityMap(bh.simmaps.b, plot=F, res=res)
bh.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ bh-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$b^h$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' 
#' 
#' ### Reconstruction of labiodental reflexes for *$p$
#' 
#' With phylogenies estimated by Chang et al.:
#+ p.simmap-c, cache=T, eval=run, results='hide'
correspondeme='p'
p.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
p.density.c <- densityMap(p.simmaps.c, plot=F, res=res)
p.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ p-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$p$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ p.simmap-b, cache=T, eval=run, results='hide'
p.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
p.density.b <- densityMap(p.simmaps.b, plot=F, res=res)
p.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ p-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$p$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' 
#' 
#' ### Reconstruction of labiodental reflexes for *$w$
#' 
#' With phylogenies estimated by Chang et al.:
#+ w.simmap-c, cache=T, eval=run, results='hide'
correspondeme='w'
w.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
w.density.c <- densityMap(w.simmaps.c, plot=F, res=res)
w.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ w-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$w$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ w.simmap-b, cache=T, eval=run, results='hide'
w.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
w.density.b <- densityMap(w.simmaps.b, plot=F, res=res)
w.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ w-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$w$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$m$
#' 
#' With phylogenies estimated by Chang et al.:
#+ m.simmap-c, cache=T, eval=run, results='hide'
correspondeme='m'
m.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ARD')
m.density.c <- densityMap(m.simmaps.c, plot=F, res=res)
m.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ARD')

#+ m-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$m$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ m.simmap-b, cache=T, eval=run, results='hide'
m.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ARD')
m.density.b <- densityMap(m.simmaps.b, plot=F, res=res)
m.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ARD')

#' 
#+ m-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$m$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$g^w$
#' 
#' With phylogenies estimated by Chang et al.:
#+ gw.simmap-c, cache=T, eval=run, results='hide'
correspondeme='gw'
gw.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
gw.density.c <- densityMap(gw.simmaps.c, plot=F, res=res)
gw.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ gw-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$g^w$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ gw.simmap-b, cache=T, eval=run, results='hide'
gw.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
gw.density.b <- densityMap(gw.simmaps.b, plot=F, res=res)
gw.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ gw-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$g^w$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$g^{wh}$
#' 
#' With phylogenies estimated by Chang et al.:
#+ gwh.simmap-c, cache=T, eval=run, results='hide'
correspondeme='gwh'
gwh.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
gwh.density.c <- densityMap(gwh.simmaps.c, plot=F, res=res)
gwh.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ gwh-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$g^{wh}$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ gwh.simmap-b, cache=T, eval=run, results='hide'
gwh.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
gwh.density.b <- densityMap(gwh.simmaps.b, plot=F, res=res)
gwh.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ gwh-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$g^{wh}$, assuming the phylogenies estimated by Bouckaert et al. ', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$k^w$
#' 
#' With phylogenies estimated by Chang et al.:
#+ kw.simmap-c, cache=T, eval=run, results='hide'
correspondeme='kw'
kw.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
kw.density.c <- densityMap(kw.simmaps.c, plot=F, res=res)
kw.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ kw-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$k^w$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ kw.simmap-b, cache=T, eval=run, results='hide'
kw.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
kw.density.b <- densityMap(kw.simmaps.b, plot=F, res=res)
kw.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ kw-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$k^w$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#+ stepcheck3, include=F
run=T
#' ### Reconstruction of labiodental reflexes for *$d^h$
#' 
#' With phylogenies estimated by Chang et al.:
#+ dh.simmap-c, cache=T, eval=run, results='hide'
correspondeme='dh'
dh.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
dh.density.c <- densityMap(dh.simmaps.c, plot=F, res=res)
dh.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ dh-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$d^h$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ dh.simmap-b, cache=T, eval=run, results='hide'
dh.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
dh.density.b <- densityMap(dh.simmaps.b, plot=F, res=res)
dh.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ dh-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$d^h$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$k^j$
#' 
#' With phylogenies estimated by Chang et al.:
#+ k.simmap-c, cache=T, eval=run, results='hide'
correspondeme='ky'
ky.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ER')
ky.density.c <- densityMap(ky.simmaps.c, plot=F, res=res)
ky.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ER')

#+ k-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$k^j$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ k.simmap-b, cache=T, eval=run, results='hide'
ky.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ER')
ky.density.b <- densityMap(ky.simmaps.b, plot=F, res=res)
ky.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ER')

#' 
#+ k-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$k^j$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' ### Reconstruction of labiodental reflexes for *$b$
#' 
#' With phylogenies estimated by Chang et al.:
#+ b.simmap-c, cache=T, eval=run, results='hide'
correspondeme='b'
b.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.sum.tree, model='ARD')
b.density.c <- densityMap(b.simmaps.c, plot=F, res=res)
b.multi.simmaps.c <- compute.simmap(correspondeme, tree=ie.c.trees.sample, model='ARD')

#+ b-scm-plot-c, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$b$, assuming the phylogenies estimated by Chang et al.', echo=F, eval=run, fig.height=11, cache=T
plot.density(correspondeme, tree=ie.c.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' With phylogenies estimated by Bouckaert et al.:
#+ b.simmap-b, cache=T, eval=run, results='hide'
b.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.sum.tree, model='ARD')
b.density.b <- densityMap(b.simmaps.b, plot=F, res=res)
b.multi.simmaps.b <- compute.simmap(correspondeme, tree=ie.b.trees.sample, model='ARD')

#' 
#+ b-scm-plot-b, fig.cap='Posterior probability (PP) density of labiodental reflexes in some contexts of the correspondeme *$b$, assuming the phylogenies estimated by Bouckaert et al.', echo=F, eval=run, fig.height=12, cache=T
plot.density(correspondeme, tree=ie.b.sum.tree, nodes=T, multi.ace=T)
#' 
#' \clearpage
#' 
#' 
#' Summary
#' -------
#' 
#' The estimated root probabilities are as follows:
#' 
#+ scm.aces-c, eval=run, echo=T, cache=T

# doing this with lapply exceeds memory, so....; and let's do it externally...
# bh.c <- as.data.frame(describe.simmap(bh.multi.simmaps.c)$ace)
# p.c <- as.data.frame(describe.simmap(p.multi.simmaps.c)$ace)
# w.c <- as.data.frame(describe.simmap(w.multi.simmaps.c)$ace)
# m.c <- as.data.frame(describe.simmap(m.multi.simmaps.c)$ace)
# gw.c <- as.data.frame(describe.simmap(gw.multi.simmaps.c)$ace)
# gwh.c <- as.data.frame(describe.simmap(gwh.multi.simmaps.c)$ace)
# kw.c <- as.data.frame(describe.simmap(kw.multi.simmaps.c)$ace)
# dh.c <- as.data.frame(describe.simmap(dh.multi.simmaps.c)$ace)
# ky.c <- as.data.frame(describe.simmap(ky.multi.simmaps.c)$ace)
# b.c <- as.data.frame(describe.simmap(b.multi.simmaps.c)$ace)

load("all_aces_c.rda")

all.aces.c <- list(bh.c, p.c, w.c, m.c, gw.c, gwh.c, kw.c, dh.c, ky.c, b.c)
names(all.aces.c) <- correspondemes					

(root.aces.c <- sapply(all.aces.c, '[', 1, 2))

mean(root.aces.c)

median(root.aces.c)

#+ scm.aces-b, eval=run, echo=T, cache=T

# bh.b <- as.data.frame(describe.simmap(bh.multi.simmaps.b)$ace)
# p.b <- as.data.frame(describe.simmap(p.multi.simmaps.b)$ace)
# w.b <- as.data.frame(describe.simmap(w.multi.simmaps.b)$ace)
# m.b <- as.data.frame(describe.simmap(m.multi.simmaps.b)$ace)
# gw.b <- as.data.frame(describe.simmap(gw.multi.simmaps.b)$ace)
# gwh.b <- as.data.frame(describe.simmap(gwh.multi.simmaps.b)$ace)
# kw.b <- as.data.frame(describe.simmap(kw.multi.simmaps.b)$ace)
# dh.b <- as.data.frame(describe.simmap(dh.multi.simmaps.b)$ace)
# ky.b <- as.data.frame(describe.simmap(ky.multi.simmaps.b)$ace)
# b.b <- as.data.frame(describe.simmap(b.multi.simmaps.b)$ace)

load("all_aces_b.rda")

all.aces.b <- list(bh.b, p.b, w.b, m.b, gw.b, gwh.b, kw.b, dh.b, ky.b, b.b)
names(all.aces.b) <- correspondemes					
					
(root.aces.b <- sapply(all.aces.b, '[', 1, 2))

mean(root.aces.b)

median(root.aces.b)


#' Figures S4.27 and S4.28 give a synopsis of all stochastic character maps. 

#+ all-maps-a, fig.cap='Stochastic character maps on the MCC trees estimated by Chang et al.', fig.height=12, eval=run, cache=T, echo=F, fig.pos='h', out.extra=''

correspondemes <- c('bh','p','w','m','gw','gwh','kw','dh','ky','b')
correspondemes.for.print <- paste0('*', c("bʰ","p","w","m","gʷ","gʷʰ","kʷ","dʰ","kʲ","b"))
par(mfrow=c(2,5), cex=.5)
x <- lapply(seq_along(correspondemes), function(i) {
    plot.density(correspondemes[i], tree=ie.c.sum.tree, nodes=F, multi.ace=F, legend=F, lwd=c(4,4))
    text(x=2, y= 28, correspondemes.for.print[i], cex=2)
    })
#' 
#' 
#+ all-maps-b, fig.cap='Stochastic character maps on the MCC trees estimated by Bouckaert et al.', fig.height=12, eval=run, echo=F, cache=T
par(mfrow=c(2,5), cex=.5)
x <- lapply(seq_along(correspondemes), function(i) {
    plot.density(correspondemes[i], tree=ie.b.sum.tree, nodes=F, multi.ace=F, legend=F, lwd=c(4,4))
    text(x=2, y= 28, correspondemes.for.print[i], cex=2)
    })

#+ include=F, echo=F
run=T
#' 
#' Figures S4.29 and S4.30 give a synopsis in terms of traitgrams. Each color represents the evolution of a correspondeme over all branches of the phylogeny. The position of the branch on the y-axis reflects the probability of a labiodental articulation.
#' 
#+ phenogram-c-prep, eval=run, cache=T
all.maps.c <- list(
    bh=bh.density.c,p=p.density.c,w=w.density.c,m=m.density.c,
    gw=gw.density.c,gwh=gwh.density.c,kw=kw.density.c,
    dh=dh.density.c,ky=ky.density.c,b=b.density.c)

all.maps.values.c <- lapply(all.maps.c, function(density.map) {
    # get node and tip probabilities from densitymap:
    nodes.and.tip.values <- c(as.numeric(names(density.map$tree$maps[[1]])[1]), # the root
            sapply(density.map$tree$maps, function(m) {
               as.numeric(names(m)[length(m)])} ) # the end states, i.e the internal nodes and tips
        )/1000
    names(nodes.and.tip.values) <- as.character(c(53, # the root number
                        density.map$tree$edge[,2]
                        )) # the internal node and tip numbers
    density.map$tree$tip.label <- as.character(1:length(density.map$tree$tip.label)) # name match
    return(list(values=nodes.and.tip.values,tree=as.phylo(density.map$tree)))
    })
    

tree.c <-  all.maps.values.c[[1]]$tree

#+ phenogram-c, fig.cap='Evolution of labiodentals across correspondemes, based on the maximum clade credibility tree of the phylogeny in Chang et al.', eval=run, cache=T,echo=F
# par(oma = c(0, 0, 1.7, 0))
x <- phenogram(tree.c, all.maps.values.c[[1]]$values,
        colors="black", 
        ftype='off', 
        ylim=c(0,1), ylab='Probability of labiodental',
        xlab='time (ybp)',
        x.axis.intervals=c(seq(0,6000,1000), round(max(nodeHeights(tree.c)),-2))
        )
print.colors <- brewer.pal(9,'Set1')
names(print.colors) <- names(all.maps.values.c)[-1]
for(i in names(all.maps.values.c)[-1]) {
    phenogram(tree.c, all.maps.values.c[[i]]$values,
        colors=print.colors[i],
        ftype='off',
        ylim=c(0,1),
        xaxt='n',
        add=T)
}
legend('top', bty='n', lwd=3, xpd = T, ncol=5, title='',  inset=c(0,-.15),
    legend = correspondemes.for.print,
    col = c('#000000',print.colors))

#' 
#' 
#+ phenogram-b-prep, eval=run, cache=T
all.maps.b <- list(
    bh=bh.density.b,p=p.density.b,w=w.density.b,m=m.density.b,
    gw=gw.density.b,gwh=gwh.density.b,kw=kw.density.b,
    dh=dh.density.b,ky=ky.density.b,b=b.density.b)

all.maps.values.b <- lapply(all.maps.b, function(density.map) {
    # get node and tip probabilities from densitymap:
    nodes.and.tip.values <- c(as.numeric(names(density.map$tree$maps[[1]])[1]), # the root
            sapply(density.map$tree$maps, function(m) {
               as.numeric(names(m)[length(m)])} ) # the end states, i.e the internal nodes and tips
        )/1000
    names(nodes.and.tip.values) <- as.character(c(53, # the root number
                        density.map$tree$edge[,2]
                        )) # the internal node and tip numbers
    density.map$tree$tip.label <- as.character(1:length(density.map$tree$tip.label)) # name match
    return(list(values=nodes.and.tip.values,tree=as.phylo(density.map$tree)))
    })
    

tree.b <-  all.maps.values.b[[1]]$tree

# /*
# plotTree(tree); nodelabels()
# plotTree(bh.density.b$tree); nodelabels()
# */

#+ phenogram-b, fig.cap='Evolution of labiodentals across correspondemes, based on the maximum clade credibility tree of the phylogeny in Bouckaert et al.', eval=run, cache=T, echo=F
# quartz(width=8.5, height=6.4)
# par(oma = c(0, 0, 1.7, 0))
# par(cex=1.3, cex.lab = 1.3)
x <- phenogram(tree.b, all.maps.values.b[[1]]$values,
        colors="black", 
        ftype='off', 
        ylim=c(0,1), 
        ylab='', # ylab='Probability of labiodental',
        xlab='', #'Time (ybp)',
        x.axis.intervals=c(seq(0,7000,2000),round(max(nodeHeights(tree.b)),-2))
        )   
print.colors <- alpha(brewer.pal(9,'Set1'),.6)
names(print.colors) <- names(all.maps.values.b)[-1]
for(i in names(all.maps.values.b)[-1]) {
    phenogram(tree.b, all.maps.values.b[[i]]$values, 
        colors=print.colors[i], 
        ftype='off', 
        ylim=c(0,1),
        xaxt='n', 
        add=T) 
}
legend('top', bty='n', lwd=3, xpd = T, ncol=5, title='',  inset=c(0,-.30), # -.15 for rmarkdown
    legend = correspondemes.for.print,
    col = c('#000000',print.colors))
title(ylab="Probability of labiodental", line=2, cex.lab=1.2)     
title(xlab="Time (ybp)", line=2.1, cex.lab=1.2)     
#' 
#' \clearpage
#' 
#' 
#' 
#' References
#' ==========
#' 