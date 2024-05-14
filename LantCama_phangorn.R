# genotype_matrix <- dms$gt
genotype_matrix <- read.nexus.data('/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen/raw_SNPFilt_1SNPperClone/svdq/LantCama_DLan23-8067_eilish.nex')
x <- nexus2DNAbin(genotype_matrix)
x2 <- as.phyDat(x)

## Not run: 
# data(x2)
dm <- dist.hamming(x2)
tree <- upgma(dm)
# NJ# NJupgma()
set.seed(123)
UPGMAtrees <- bootstrap.phyDat(x2,
                            FUN=function(x)upgma(dist.hamming(x)), bs=100)
treeUPGMA <- plotBS(tree, UPGMAtrees, "phylogram")

# Maximum likelihood
fit <- pml(tree, x2)
fit <- optim.pml(fit, rearrangement="NNI")
set.seed(123)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)

# Maximum parsimony
treeMP <- pratchet(x2)
treeMP <- acctran(treeMP, x2)
set.seed(123)
BStrees <- bootstrap.phyDat(x2, pratchet, bs = 100)
treeMP <- plotBS(treeMP, BStrees, "phylogram")
add.scale.bar()

# export tree with bootstrap values as node labels
# write.tree(treeBS)

## End(Not run)