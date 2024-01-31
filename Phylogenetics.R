# Evolutionary Genetics and Genomics (HGEN 6092)
# R Session 3 - Phylogenetics
# Casey Meili

# Some exercises were adapted from: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

# load packages
library(ape)
library(adegenet)
library(phytools)
library(phangorn)

# import downloaded dataset
dna <- fasta2DNAbin(file="/Users/caseymeili/Downloads/PIP_primate.fa")  

# DISTANCE BASED APPROACHES
# dist.dna() from the {ape} package makes a matrix of pairwise distances from the DNA sequences. “TN93” is the evolutionary model - this one allows for different transition rates and heterogenous base frequencies
# compute genetic distance between all pairs of sequences and visualize 
D <- dist.dna(dna, model = "TN93") 
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) # darker shades of gray mean a larger distance

# build neighbor joining tree using computed distances
tre <- nj(D)
class(tre) 
tre <- ladderize(tre)
tre # output tells us what the tree will look like but doesn't show the actual construction

plot(tre,type="unrooted")
title("Unrooted NJ tree")

# rooting the tree using bushbaby as an outgroup
tre2 <- root(tre, out = "bushbaby")
tre2 <- ladderize(tre2)
plot(tre2)
title("Rooted NJ tree")
axisPhylo()  # adds a scaled axis of distances

# bootstrapping the phylogeny
myBoots <- boot.phylo(tre2, dna, function(e) root(nj(dist.dna(e, model = "TN93")),"bushbaby"))
plot(tre2, edge.width=2)
title("NJ tree + bootstrap values")
nodelabels(myBoots, cex=.6)

# MAXIMUM PARSIMONY APPROACH
# import data and assign as object
dna2 <- as.phyDat(dna) # assign the original dna sequences data as a phyDat object
dna2

# reconstruct the tree using neighbor-joining, as before
tre.ini <- nj(dist.dna(dna,model="raw")) 
tre.ini

# compute a parsimony score for the tree (~how many changes it requires)
parsimony(tre.ini, dna2)  

# make the most parsimonious tree possible
tre.pars <- optim.parsimony(tre.ini, dna2) 

# compute a parsimony score
parsimony(tre.pars, dna2)

# plot maximum parsimony tree
plot(tre.pars, type="unr", show.tip=TRUE, edge.width=2)
title("Maximum-parsimony tree")

# MAXIMUM LIKELIHOOD APPROACH
# calculate the likelihood of the data given the model for the neighbor-joining tree
fit.ini <- pml(tre.ini, dna2, k=4)
fit.ini

# optimize the tree
fit <- optim.pml(fit.ini, optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)
fit

# compare neighbor joining to maximum likelihood tree using anova
anova(fit.ini, fit)

# compare trees using akaike information criterion (AIC), lower AIC value indicates better fit
AIC(fit.ini)
AIC(fit)

# plot ML tree and root using the same outgroup
tre4 <- root(fit$tree,"bushbaby")
tre4 <- ladderize(tre4)
plot(tre4, edge.width=2)
title("Maximum-likelihood tree")  

