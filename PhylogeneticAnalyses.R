# N.B. Do not run the computational intense analyses in RStudio (e.g., RAxML, MrBayes)

#########################
# Define useful variables
#########################
directory <- '~/Desktop/MolecularEvolutionAnalyses/'

#######################
# Set working directory
#######################
setwd(directory)

##################
# Install packages
##################
install.packages('ape', dependencies = TRUE)
install.packages('ips', dependencies = TRUE)
install.packages('phangorn', dependencies = TRUE)
install.packages('spatstat', dependencies = TRUE)
install.packages('expm', dependencies = TRUE)
library(ape)
library(ips)
library(phangorn)
library(spatstat)
library(expm)

###############################
# Import data (with outgroups!)
###############################
COI.all <- read.dna(file = 'COI.all.fas', format = 'fasta')
r16S.all <- read.dna(file = '16S.all.fas', format = 'fasta')
ITS2.all <- read.dna(file = 'ITS2.all.fas', format = 'fasta')

# TBD: import to three fasta files at once

## Extract labels
COI.all.labels <- labels(COI.all)
r16S.all.labels <- labels(r16S.all)
ITS2.all.labels <- labels(ITS2.all)

# TBD: extract the labels of the three fasta files at once

## Find the common sequences
intersect(intersect(COI.all.labels, r16S.all.labels), ITS2.all.labels)

###########################
# Exploratory data analyses
###########################
exploratory.data.analyses <- function(fas){
  structure <- str(fas)
  names <- names(fas)
  labels <- labels(fas)
  num.seq <- length(fas)
  list <- list(structure, names, labels, num.seq)
  return(list)
}

exploratory.data.analyses(COI.all)
exploratory.data.analyses(r16S.all)
exploratory.data.analyses(ITS2.all)

## Base frequency
base.freq(COI.all)
#base.freq(ITS2.all)[c('-', 'n')]
### Base frequency for each individual
BF.individual <- function(fas){
  BF <- matrix(NA, nrow = length(fas), ncol = 4)
  rownames(BF) <- fas$labels
  colnames(BF) <- c('A', 'C', 'G', 'T')
  for (i in 1:length(fas)) {
    BF[i, ] <- base.freq(fas[i])
  }
  plot <- matplot(BF, type = 'l', col = 1, xlab = 'Species', ylab = 'Base frequency')
  legend(x = 'left', c('A', 'C', 'G', 'T'), lty=1:4, bty = 'n')
}

BF.individual(COI.all)
BF.individual(r16S.all)
BF.individual(ITS2.all) # TBD: better position the legend

### Base frequency between genes for all species pooled (p. 74) - WIP
common.individuals <- intersect(intersect(COI.all.labels, r16S.all.labels), ITS2.all.labels)
COI.com <- COI.all[common.individuals]
r16S.com <- r16S.all[common.individuals]
ITS2.com <- ITS2.all[common.individuals]

## Site rate heterogeneity - WIP

#############################
# Multiple sequence alignment
#############################

## CLUSTAL
COI.all.clustal <- clustal(COI.all)
r16S.all.clustal <- clustal(r16S.all)
ITS2.all.clustal <- clustal(ITS2.all)

# TBD: align the three fasta files at once

### Explore the alignment
explore.alignment <- function(alignment){
  dimensions <- dim(alignment)
  image(alignment)
  image(alignment, 'a') 
  image(alignment, '-')
  image(alignment, c('g', 'n'), c('black', 'grey'))
  grid(ncol(alignment), nrow(alignment), col = 'lightgrey')
  checkAlignment(alignment)
  alview(alignment)
  return(dimensions)
}

explore.alignment(COI.all.clustal)
explore.alignment(r16S.all.clustal)
explore.alignment(ITS2.all.clustal)

# How is the Shannon index for nucleotide diversity defined?

## MUSCLE
r16S.all.muscle <- muscle(r16S.all)
ITS2.all.muscle <- muscle(ITS2.all)

# TBD: align the two fasta files at once

### Explore the alignment
explore.alignment(r16S.all.muscle)
explore.alignment(ITS2.all.muscle)

## MAFFT
r16S.all.mafft <- mafft(r16S.all)
ITS2.all.mafft <- mafft(ITS2.all)

# TBD: align the two fasta files at once

### Explore the alignment
explore.alignment(r16S.all.mafft)
explore.alignment(ITS2.all.mafft)

dev.off(dev.list()['RStudioGD'])

## Compare the alignments
all.equal.DNAbin(r16S.all.clustal, r16S.all.muscle, plot = TRUE)
all.equal.DNAbin(r16S.all.clustal, r16S.all.mafft, plot = TRUE)
all.equal.DNAbin(r16S.all.muscle, r16S.all.mafft, plot = TRUE)
# Why is there a gap in position 18 and not 19?

all.equal.DNAbin(ITS2.all.clustal, ITS2.all.muscle, plot = TRUE)
all.equal.DNAbin(ITS2.all.clustal, ITS2.all.mafft, plot = TRUE)
all.equal.DNAbin(ITS2.all.muscle, ITS2.all.mafft, plot = TRUE)

dev.off(dev.list()['RStudioGD'])

## Cut the sequences where the forward primer ends and where the reverse primer starts - WIP
COI.all.clustal <- COI.all.clustal[ , 1:710]
r16S.all.clustal <- r16S.all.clustal[ , 2:547]
r16S.all.muscle <- r16S.all.muscle[ , 2:548]
r16S.all.mafft <- r16S.all.mafft[ , 2:548]
ITS2.all.clustal <- ITS2.all.clustal[ , 2:351]
ITS2.all.muscle <- ITS2.all.muscle[ , 2:354]
ITS2.all.mafft <- ITS2.all.mafft[ , 2:354]

## Discard duplicate haplotypes
#COI.all.clustal <- as.phyDat(COI.all.clustal)
#COI.all.clustal <- unique(COI.all.clustal) # Alternative : haplotype() from ape
#COI.all.clustal <- as.DNAbin(COI.all.clustal)

## Write alignmenents
#write.dna(COI.all.clustal, 'COI.all.clustal.fas', 'fasta', colsep = '')
#write.dna(r16S.all.mafft, 'r16S.all.clustal.fas', 'fasta', colsep = '')
#write.dna(r16S.all.mafft, 'r16S.all.muscle.fas', 'fasta', colsep = '')
#write.dna(r16S.all.mafft, 'r16S.all.mafft.fas', 'fasta', colsep = '')
#write.dna(ITS2.all.mafft, 'ITS2.all.mafft.fas', 'fasta', colsep = '')

# Replace external gaps by Ns in AliView

############################################################################################################################

## Import final alignments
COI.all.clustal <- read.dna(file = 'COI.all.clustal.fas', format = 'fasta')
r16S.all.mafft <- read.dna(file = 'r16S.all.mafft.fas', format = 'fasta')
# In my research group, we usually use the secondary structure to align 16S. 
# I thought I would also align it with algorithms to check the differences with the secondary structure alignment.
# Why did I choose MAFFT?
ITS2.all.mafft <- read.dna(file = 'ITS2.all.mafft.fas', format = 'fasta') 
# I chose the MAFFT alignment for ITS2 because it is the algorithm which is used in my research group for this gene and for gastropods.
# But why?  

### Extract labels
COI.all.clustal.labels <- labels(COI.all.clustal)
r16S.all.mafft.labels <- labels(r16S.all.mafft)
ITS2.all.mafft.labels <- labels(ITS2.all.mafft)

# TBD: extract the labels of the three fasta files at once

### Find common sequences
common.individuals.mt <- intersect(COI.all.clustal.labels, r16S.all.mafft.labels)
COI.mt <- as.DNAbin(subset(as.phyDat(COI.all.clustal), common.individuals.mt))
r16S.mt <- as.DNAbin(subset(as.phyDat(r16S.all.mafft), common.individuals.mt))

#write.dna(COI.mt, 'COI.fas', 'fasta', colsep = '')
#write.dna(r16S.mt, '16S.fas', 'fasta', colsep = '')

identical(labels(COI.mt), labels(r16S.mt))

## The class 'DNAbin'
COI.all.clustal[1:5, 1:100]
COI.all.clustal[sort(rownames(COI.all.clustal)), ]
# To look at the third position of each codon, position which evolves faster (?) than the two other ones. 
# Only makes sense for protein-coding sequences. 
s <- c(FALSE, FALSE, TRUE) 
COI.all.clustal[ , s]
COI.all.clustal[ , !s] # ! inverts the logical value
image(COI.all.clustal[ , s])

dev.off(dev.list()['RStudioGD'])

## Calculate the base frequency between sites for a single gene (p. 76)
table(sapply(COI.all.clustal, length)) # looking at the length of each sequence
nm <- rownames(COI.all.clustal)
BF.COI.all.clustal <- matrix(NA, 3, 4)
rownames(BF.COI.all.clustal) <- paste('codon position', 1:3)
colnames(BF.COI.all.clustal) <- c('A', 'C', 'G', 'T')
for (i in 1:3) {
  s <- rep(FALSE, 3)
  s[i] <- TRUE
  BF.COI.all.clustal[i, ] <- base.freq(COI.all.clustal[, s])
}
BF.COI.all.clustal
barplot(t(BF.COI.all.clustal), main = 'Cytochrome oxydase I', ylab = 'Base frequency')
par(cex = 2)
text(0.7, BF.COI.all.clustal[1, 1]/2, 'A', col = 'white')
text(0.7, BF.COI.all.clustal[1,1]+BF.COI.all.clustal[1,2]/2, 'C', col = 'white')
text(0.7, sum(BF.COI.all.clustal[1, 1:2]) + BF.COI.all.clustal[1, 3]/2, 'G')
text(0.7, sum(BF.COI.all.clustal[1, 1:3]) + BF.COI.all.clustal[1, 4]/2, 'T')

dev.off(dev.list()['RStudioGD'])

## Homogeneity partition test - WIP

## Relative apparent synapomorphy analysis + permutation probability test - WIP

## Power test - WIP

#########################################
# Best-fit substitution models estimation
#########################################
model.test <- function(alignment){
  alignment.phyDat <- as.phyDat(alignment)
  mt <- modelTest(alignment.phyDat, multicore = TRUE, mc.cores = 4) # Using multicore processors does not work on Windows
  # phymltest {ape} is no longer maintained
  best.model.AIC <- mt$Model[which.min(mt$AIC)] # Choosing the best model according to AIC
  best.model.AICc <- mt$Model[which.min(mt$AICc)] # Choosing the best model according to AIC
  best.model.BIC <- mt$Model[which.min(mt$BIC)] # Choosing the best model according to AIC
  if (is.empty(best.model.AICc)) {
    df <- data.frame(Criterion = c('AIC', 'BIC'),
                     Model = c(best.model.AIC, best.model.BIC))
  } else {
    df <- data.frame(Criterion = c('AIC', 'AICc', 'BIC'),
                     Model = c(best.model.AIC, best.model.AICc, best.model.BIC))
  }
  return(df)
}

model.test.COI <- model.test(COI.all.clustal)
model.test.COI
model.test.r16S <- model.test(r16S.all.mafft)
model.test.r16S
model.test.ITS2 <- model.test(ITS2.all.mafft)
model.test.ITS2

#################
# Saturation test
#################
saturation <- function(alignment){
  TS <- dist.dna(alignment, 'ts')
  TV <- dist.dna(alignment, 'tv')
  tree <- dist.dna(alignment, 'K80') 
  if (max(TV) < max(TS)){
    plot(tree, TS, ylim = c(0, max(TS)))
  }else{
    plot(tree, TS, ylim = range(TV, TS))
  }
  points(tree, TV, pch = 2)
  legend('topleft', c('Transition', 'Transversion'), pch = 1:2)
}

saturation(COI.all.clustal)
saturation(r16S.all.mafft)
saturation(ITS2.all.mafft)
# Interpretation?

# See Wilke et al. 2011 for saturation test in Hydrobiidae

############################################################################################################################

## Import final alignments
# I import them again because the analyse below were not working with the previously imported alignments? 
# But why? Both series of alignments were of the class 'DNAbin'
COI.all.clustal <- read.FASTA('COI.all.clustal.fas')
r16S.all.mafft <- read.FASTA('r16S.all.mafft.fas')
ITS2.all.mafft <- read.FASTA('ITS2.all.mafft.fas') 

### Extract labels
COI.all.clustal.labels <- labels(COI.all.clustal)
r16S.all.mafft.labels <- labels(r16S.all.mafft)
ITS2.all.mafft.labels <- labels(ITS2.all.mafft)

# TBD: extract the labels of the three fasta files at once

### Find common sequences
mt <- intersect(COI.all.clustal.labels, r16S.all.mafft.labels)
COI.mt <- COI.all.clustal[mt]
r16S.mt <- r16S.all.mafft[mt]
identical(labels(r16S.mt), labels(COI.mt))

############################
# Merge the alignments - WIP
############################

mt.genes <- list()
n <- length(labels(COI.mt))
for(i in 1:n){
  mt.genes[[i]] <- c(COI.mt[[i]], r16S.mt[[i]])
}
class(mt.genes) <- 'DNAbin'
image(mt.genes)
class(mt.genes)
str(mt.genes)
labels(mt.genes)
names(mt.genes) <- labels(COI.mt)
labels(mt.genes)

#############
# Infer trees
#############

### Brute force
howmanytrees(length(mt.genes.labels), rooted = FALSE)

## Evolutionary distances
Q <- matrix(c(-0.1, 0.2, 0.1, -0.2), 2) # Rates estimated based on data
matexpo(Q)
matexpo(Q*100)
expm(Q*100)

## Distance-based phylogenetic methods
draw <- dist.dna(mt.genes, 'raw') # 'Raw' distances
djc69 <- dist.dna(mt.genes, 'JC69')
dk81 <- dist.dna(mt.genes, 'K81')
plot(djc69, dk81)
# Interpretation?

### Possible diagnostic to help deciding whether or not using pairwise deletion
# What are pairwise deletion?
sum(base.freq(mt.genes, freq = TRUE, all = TRUE)[c('-', 'n')]) # If = 0, no need to use pairwise.deletion = TRUE

#### Compare pairwise deletion
summary.default(dk81 - dist.dna(mt.genes, pairwise.deletion = TRUE))
summary.default(djc69 - dist.dna(mt.genes, pairwise.deletion = TRUE))
plot(dk81, dist.dna(mt.genes, pairwise.deletion = TRUE))
abline(0, 1)
# Interpretation?

### Check the distribution of distances
layout(matrix(1:2, 1))
hist(dist.dna(mt.genes), main = 'dist'); rug(dist.dna(mt.genes))
hist(dk81, breaks = 50, main = 'dk81'); rug(dk81)
# Interpretation?

### Assess the treelikeliness of phylogenetic distance data before tree estimation (see Holland et al. 2002)
delta.plot(djc69)
delta.plot(dk81)
all(attr(djc69, 'Labels') == attr(dk81, 'Labels'))
# Interpretation?

dev.off(dev.list()['RStudioGD'])

### Number of polymorphic/segregating sites
seg.sites(mt.genes)
#Segregating sites are positions which show differences (polymorphisms) between related genes in a sequence alignment (are not conserved).
#Segregating sites include conservative, semi-conservative and non-conservative mutations.
#The proportion of segregating sites within a gene is an important statistic in population genetics since it can be used to estimate mutation rate assuming no selection. 
#For example it is used to calculate the Tajima's D neutral evolution statistic (Wikipedia).

### UPGMA
tr.upgma <- upgma(djc69)
plot(tr.upgma)

### Neighbor-Joining
dist.dna.mt.genes <- dist.dna(mt.genes)
tr.nj <- nj(dist.dna.mt.genes)
plot(tr.nj)
tr.nj.rooted <- root(tr.nj, c('Safe2231', 'Peul1968'))
# tr.nj.rooted <- drop.tip(tr.nj.rooted, c('Safe2231', 'Peul1968'))
plot(tr.nj.rooted)

#### Plot the difference between distances from the estimated tree and original distances (distance residuals) for four tree estimation methods
dt.nj <- cophenetic(tr.nj)
dmat <- as.matrix(mt.genes)
nms <- rownames(dmat)
dt.nj <- dt.nj[nms, nms]
dt.nj <- as.dist(dt.nj)
plot(dt.nj - dist.dna.mt.genes, ylab = 'distance residuals')
abline(h = 0, lty = 3)
# Interpretattion? 

### Bootstrap
args(boot.phylo)
BP <- boot.phylo(tr.nj, as.DNAbin(as.phyDat(mt.genes)), function(x) nj(dist.dna(x))) 
plot(tr.nj, 'p', FALSE) # Simple plot
drawSupportOnEdges(BP)
#### Attach the bootstrap values to the tree (as node label)
tr.nj$node.label <- BP 
#### Root the tree and use the option to keep the labels attached to the edges (and not to the nodes)
tr.nj.boot <- root(tr.nj, c('Safe2231', 'Peul1968'), edgelabel = TRUE) # Does not make sense to bootstrap a rooted tree
plot(tr.nj.boot, show.node.label = TRUE) # Simple plot with bootstrap values
plot(tr.nj.boot)
nodelabels(tr.nj.boot$node.label, adj = 1)
#### Repeat the bootstrap and get the bootstrap trees
BP <- boot.phylo(tr.nj, as.DNAbin(as.phyDat(mt.genes)), function(x) nj(dist.dna(x)), 1000, trees = TRUE)
lento(BP$trees[1:1000]) # Lento plot
lento(BP$trees[1:100])
# Why are the lento plots almost completely black? 
CN <- consensusNet(BP$trees[1:1000])
# What is the difference between a consensus phylogentic tree and an haplotype network?
plot(CN)
# Why is there a fatal error when I quit Xquartz after plotting the consensus net? 

### Minimum Evolution
tr.me.bal <- fastme.bal(dist.dna.mt.genes)
plot(tr.me.bal)
tr.me.bal.rooted <- root(tr.me.bal, c('Safe2231', 'Peul1968'))
plot(tr.me.bal.rooted)
tr.me.ols <- fastme.ols(dist.dna.mt.genes)
plot(tr.me.ols)

### Check whether or not trees are the same
all.equal(tr.nj, tr.me.bal)
dist.topo(tr.nj, tr.me.bal)
# Interpretation? 
dist.topo(tr.nj.rooted, tr.me.bal.rooted)
layout(matrix(1:2, 1))
plot(tr.nj, main = 'NJ')
#add.scale.bar()
plot(tr.me.bal, main = 'ME')
#add.scale.bar()

dev.off(dev.list()['RStudioGD'])

## Maximum parsimony

### Tree rearrangements
x <- as.phyDat(mt.genes)
tree <- nj(dist.logDet(x))
treeNNI <- optim.parsimony(tree, x, rearrangements = 'NNI', trace = 0)
treeSPR <- optim.parsimony(tree, x, trace = 0)
treeRatchet <- pratchet(x, trace = 0)
fitch(c(treeNNI, treeSPR, treeRatchet), x)
as.matrix(RF.dist(c(treeNNI, treeSPR, treeRatchet)))
treeNNI <- acctran(treeNNI, x)
plot(treeNNI)
attr(x, 'index')[15]

### Ancestral reconstruction
anc.p <- ancestral.pars(treeRatchet, x)

## Maximum Likelihood methods

### Estimation with molecular sequences
dist.mt.genes <- dist.dna(mt.genes)
tree <- nj(dist.mt.genes)
data <- as.phyDat(mt.genes)
bf <- base.freq(mt.genes)
methods(class = 'pml')
#update() # Alternative
tr.ml.HKY.G.I <- pml(tree = tree, data = data, bf = bf, model = 'HKY + G + I', optInv = TRUE, optGamma = TRUE) 
tr.ml.HKY.G.I

### Optimize the model
control <- pml.control(trace = 1)
tr.ml.HKY.G.I <- optim.pml(tr.ml.HKY.G.I, Q = c(1, 1, 1, 1), model = 'HKY', optInv = TRUE, optGamma = TRUE, control = control, rearrangement = 'stochastic', subs = c(1, 1, 1, 1))
tr.ml.HKY.G.I
tr.ml.HKY.G.I$rate

### Comparing several models
tr.ml.HKY.I <- pml(tree = tree, data = data, bf = bf, model = 'HKY + I', optInv = TRUE) 
tr.ml.HKY.I
tr.ml.HKY.I <- optim.pml(tr.ml.HKY.I, Q = c(1, 1, 1, 1), model = 'HKY', optInv = TRUE, control = control, rearrangement = 'stochastic', subs = c(1, 1, 1, 1))
tr.ml.HKY.I
tr.ml.HKY.I$rate

anova(tr.ml.HKY.G.I, tr.ml.HKY.I)  # Comparing the different models

tr.ml.HKY.I.rooted <- root(tr.ml.HKY.I$tree, c('Safe2231', 'Peul1968'))
plot(tr.ml.HKY.I.rooted)
?trex

### Bootstrap analysis & visualization for ml trees
bs <- bootstrap.pml(tr.ml.HKY.I, bs = 1000, optNni = TRUE)
plotBS(tr.ml.HKY.I$tree, bs, type = 'phylo', bs.adj = c(1.5, 0), cex = 0.6)
cnet <- consensusNet(bs)
plot(cnet)

### Partition models (p. 154)
#### Set up partition models with a list of pml objects
#### Partition via a weight matrix
#### Read multiple alignments with the apex package

### Find the maximum likelihood tree (with RAxML ?)
exec <- '~/Desktop/MolecularEvolutionAnalyses/RAxMLR/raxmlHPC-PTHREADS-AVX' # Name of the executable and not of the path
partitions <- raxml.partitions(COI = as.DNAbin(as.phyDat(COI.mt)),
                               r16S = as.DNAbin(as.phyDat(r16S.mt)))
#raxml(as.DNAbin(as.phyDat(mt.genes)), m = 'GTRCAT', N = 100, b = 1234, outgroup = c('Safe2231', 'Peul1968'), partitions = partitions, backbone = tr.nj, file = 'fromR', exec = exec, threads = 4)

#### RAxML best-known tree
raxml.tree <- read.tree('~/Desktop/MolecularEvolutionAnalyses/RAxMLR/RAxML_bestTree.fromR')
plot(raxml.tree)
nodelabels(raxml.tree$node.label, adj = c(1, 0), frame = 'none')
#### RAxML bootstrap trees
raxml.bootstrap <- read.tree('RAxML_bootstrap.fromR')

## Bayesian inference

### Species tree from gene tree - WIP

#### Change individual labels to species labels
mt.genes.all.species <- mt.genes
mt.genes.all.species.labels <- labels(mt.genes.all.species)

for (i in 1:length(mt.genes.all.species.labels)){
  names(mt.genes.all.species)[i] <- substr(mt.genes.all.species.labels[i], 1, 4)
}

### Bayesian posterior inference - WIP

### Examine convergence

directory.MCMC <- '/Users/jvddorpe/Desktop/MolecularEvolutionAnalyses/AnanlysesRanByBjörn/Ecrobia\ 4/starBEAST '
install.packages('devtools')
library(devtools)
install_github('danlwarren/RWTY')
library(rwty)

#### Set the number of cores
rwty.processors <<- 4
#### Load trees
my.starbeast.trees <- load.trees('/Users/jvddorpe/Desktop/MolecularEvolutionAnalyses/AnanlysesRanByBjörn/Ecrobia\ 4/starBEAST/ecrobia4_rel_BD_starBEAST.species.trees', format = "*beast")
#### Analyze RWTY
ecrobia.rwty <- analyze.rwty(my.starbeast.trees, burnin = 50, fill.color = 'likelihood')
names(ecrobia.rwty) # See the plots I have
#### Make the plots
makeplot.all.params(my.starbeast.trees, burnin=0)
makeplot.all.params(my.starbeast.trees, burnin=50)

## Codon model

## The class 'phylo'

## Visualizaton methods

## Plot trees

#################
# Tree comparison
#################

###########
# Bootstrap
###########

############################
# Divergence time estimation
############################
## Likelihood ratio test
## Calculation of uncorrected average pairwise sequence divergence (dxy) for the split of lineages
## Calculation of standard error (SE) for the split lineages
## Population divergence (tau)
## Correction for within-region polymorphism
## Selecting the model of lineage specific rate variation
## Calibration
### Biogeographical
### Birth-death process

## Root

#####################
# Coalescence analyis
#####################

############
# List trees
############
## The class 'multiphylo'

############
# Tree space
############

##########
# phytools
##########