library(phangorn)
library(ape)

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
                            FUN=function(x)upgma(dist.hamming(x)), bs=1000)
treeUPGMA <- plotBS(tree, UPGMAtrees, "phylogram")

write.tree(treeUPGMA, file='treeUPGMA.tree')
write.nexus(treeUPGMA, file='treeUPGMA.nex')

save(tree, treeUPGMA, file = "/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/treeUPGMA.RData")


# Maximum likelihood
fit <- pml(tree, x2)
fit <- optim.pml(fit, rearrangements="stochastic")
set.seed(123)
bs <- bootstrap.pml(fit, bs=1000, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)
save(fit, bs, file = "/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/treeML.RData")

# Maximum parsimony
treeMP <- pratchet(x2)
treeMP <- acctran(treeMP, x2)
set.seed(123)
BStrees <- bootstrap.phyDat(x2, pratchet, bs = 1000)
treeMP <- plotBS(treeMP, BStrees, "phylogram")
add.scale.bar()
save(treeMP, BStrees, file = "/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/treeMP.RData")

# export tree with bootstrap values as node labels
# write.tree(treeBS)

## End(Not run)


########

upgma <- read.tree('treeUPGMA.tree')

ggtree_obj <- ggtree(upgma)

x1 <- m2[m2$sample %in% ggtree_obj$data$label,]
rownames(x1) <- x1$sample


ggtree_obj <- ggtree(upgma, size=0.02) %<+% x1 

cc <-   named_list_maker(x1$svdq_pop, 'Spectral',11)
cc2 <- named_list_maker(x1$national2, 'Paired',11)

hmt <-  gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop','national2')]),
                 offset=0.001, width=.05,font.size=2,
                 colnames_angle=90, colnames_position="top",
                 custom_column_labels=c("SVDq cluster","Country"), hjust=0)+
  scale_fill_manual(values=c(cc,cc2), na.value = "grey90")+
  theme_tree2()+
  geom_rootedge(0.0005, size=0.01)+
  geom_tiplab(aes(label = label), size=0.75)+
  geom_label2(data=upgma,aes(label=label,
                             subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
              color='red',nudge_x = 0.00015, label.size=0, fill="transparent", size=0.75)+
  theme(legend.position = "right")+
  scale_y_continuous(expand=c(0.05, 0))


library(ggtree)
library(ggplot2)

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop','national2')]),
                offset=0.001, width=.05, font.size=2,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("SVDq cluster","Country"), hjust=0) +
  scale_fill_manual(values=c(cc,cc2), na.value = "grey90") +
  theme_tree2() +
  geom_rootedge(0.0005, size=0.01) +
  geom_tiplab(aes(label = label), size=0.75) +
  geom_label2(data=upgma, aes(label=label,
                              subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
              color='red', nudge_x = 0.00015, label.size=0, fill="transparent", size=0.75) +
  theme(legend.position = "none") +
  scale_y_continuous(expand=c(0.05, 0))
# guides(fill = guide_legend(title = "Legend Title", 
#                            override.aes = list(fill = c("svdq_pop", "national2"))))


# hmt
ggsave("LantCama/outputs/LantCama_hamming_bootstrapped_upgma_maf2.pdf",
       hmt, width = 20, height = 30, units = "cm", dpi=600)

library(ggtree)
library(ggplot2)
library(ggtree)
library(ggplot2)

# Create separate ggplots for each fill scheme
plot_svdq_pop <- ggplot(m2, aes(x=long,y=lat, fill=svdq_pop)) +
  geom_tile() +
  scale_fill_manual(name = "SVDq cluster", values = cc, na.value = "grey90") 

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = cc2, na.value = "grey90")

# Extract legends from the ggplots
legend_svdq_pop <- cowplot::get_legend(plot_svdq_pop)
legend_national2 <- cowplot::get_legend(plot_national2)

legends <- cowplot::plot_grid(legend_svdq_pop, legend_national2, ncol=1, axis='bt') + theme(aspect.ratio = 3/1) # Adjust aspect ratio as needed
# Create ggtree plot and add extracted legends
hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop','national2')]),
                offset=0.001, width=.05, font.size=2,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("SVDq cluster","Country"), hjust=0) +
  theme_tree2() +
  geom_rootedge(0.0005, size=0.01) +
  geom_tiplab(aes(label = label), size=0.75) +
  geom_label2(data=upgma, aes(label=label,
                              subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
              color='red', nudge_x = 0.00015, label.size=0, fill="transparent", size=0.75) +
  theme(legend.position = "right") +
  scale_y_continuous(expand=c(0.05, 0))

legends <- cowplot::plot_grid(legend_svdq_pop, legend_national2, ncol=1, axis="b")

# Combine the ggtree plot with the extracted legends
combined_plot <- cowplot::plot_grid(hmt, legends,nrow = 1, rel_widths = c(1, 0.2))


ggsave("LantCama/outputs/LantCama_hamming_bootstrapped_upgma_maf2.pdf",
       combined_plot, width = 20, height = 30, units = "cm", dpi=600)
