library(patchwork) # For combining plots
library(ade4)
library(adegenet)
library(ape)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(dbscan)
library(diveRsity)
library(dplyr)
library(fastDiversity)
library(geosphere)
library(ggforce)
library(ggfortify)
library(ggmap)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsignif)
library(ggthemes)
library(ggtree)
library(heatmaply)
library(lattice)
library(magrittr)
library(multipanelfigure)
library(openxlsx)
library(ozmaps)
library(RColorBrewer)
library(RRtools)
library(Rtsne)
library(SNPRelate)
library(stringr)
library(tanggle)
library(tidyr)

source('https://github.com/eilishmcmaster/SoS_functions/blob/33bd7065baa91b7e2a1800b6481d63573fb38d28/dart2svdquartets.r?raw=TRUE')
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory 
missingness <- 0.3
species   <- "LantCama"
dataset   <- "DLan23-8067"

basedir <- ""

d1        <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)
meta      <- read.meta.data.full.analyses.df(d1, basedir, species, dataset)

d3        <- dart.meta.data.merge(d1, meta) 

EA_only <- c(meta$sample_names[which(!is.na(meta$analyses[,'EA_AM']))])

d3.1 <- remove.by.list(d3,EA_only)
d4        <- remove.poor.quality.snps(d3.1, min_repro=0.98, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
length(d5$locus_names)
dms <- remove.by.missingness(d5, 0.8)
length(dms$sample_names)

m2 <- dms$meta$analyses %>% as.data.frame
m2$lat <- as.numeric(m2$lat)
m2$long <- as.numeric(m2$long)

#### read clusters ###

m_clusters <- read.xlsx('LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx')
m2$clusters <- m_clusters$cluster[match(m2$sample, m_clusters$sample)] %>% as.vector()

#### colours 
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")


n_clusters <- length(unique(m2$cluster[!is.na(m2$cluster)]))  # Exclude NA from the count
tsne_cols <- brewer.pal(n_clusters, "Paired")
tsne_cols <- c("white", tsne_cols)
names(tsne_cols) <- c(NA, unique(m2$cluster[!is.na(m2$cluster)]))
tsne_cols2 <- tsne_cols[!is.na(names(tsne_cols))]

#### splitstree ####
dms_maf2 <- remove.by.maf(dms, 0.02)
# splitstree(dist(dms_maf2$gt, method = "euclidean"), 'LantCama/outputs/20241205_all_LantCama_nexus_file_for_R_80missing_maf2.nex')

library(tanggle)
library(RSplitsTree)
# splitstree(dist(dms$gt, method = "euclidean"), 'LantCama/outputs/20241205_all_LantCama_nexus_file_for_R_80missing.nex')
library(ggplot2)
library(phangorn)
library(ggforce)

# Read network data from Nexus file
Nnet <- phangorn::read.nexus.networx('LantCama/outputs/20241205_all_LantCama_nexus_file_for_R_80missing_maf2.nex')


x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)
x$cluster[!is.na(x$sample)&is.na(x$cluster)] <- "ungrouped"
# 
net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)


levels2 <- unique(dms$meta$analyses[,"national2"]) %>% sort()
levels2_shapes <- setNames(1:nlevels(factor(levels2)), levels2)

hull <- x %>% group_by(cluster) %>%
  slice(chull(x, y))
hull <- hull[!hull$cluster=='ungrouped',]

morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")

#


# Create hull data
hull <- x %>%
  group_by(cluster) %>%
  slice(chull(x, y)) %>%
  filter(!cluster == "ungrouped")  # Exclude ungrouped clusters

# Filter points for labeling
labels <- hull %>%
  group_by(cluster) %>%
  filter(
    (y > 20 & y == max(y)) |
           # (y > 20  & x > -50 & x == min(x)) |
           (y < -10 & x < -20 & y == min(y) )|
           (y > -10 & x > -20 & x == max(x) )
    ) %>%
  mutate(y = if_else(y > -10, y + 1, y - 5)) %>%
  # mutate(x = if_else((y > -10 & x < -50), x + 2, y + 0)) %>%
  ungroup()

# Add labels to the plot
splitstree_plot <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_splitnet(layout = "slanted", size = 0.2) +
  scale_fill_manual(values = tsne_cols, na.translate = FALSE) +
  scale_colour_manual(values = morphid_colours, na.translate = FALSE) +
  theme_void() +
  expand_limits(x = c(min(x_filtered$x) - 0.01 * net_x_axis, max(x_filtered$x) + 0.1 * net_x_axis),
                y = c(min(x_filtered$y) - 0.01 * net_y_axis, max(x_filtered$y) + 0.01 * net_y_axis)) +
  theme(legend.position = "bottom", legend.key = element_blank(),
        legend.key.size = unit(0, 'lines')) +
  coord_fixed() +
  labs(color = "Morphotype", shape = "Origin", fill = "HDBSCAN clusters") +
  guides(colour = guide_legend(title.position = "top", nrow = 2, override.aes = list(fill = NA, linetype = 0)), 
         fill = guide_legend(title.position = "top", nrow = 2, override.aes = list(color = NA)),
         shape = guide_legend(title.position = "top", nrow = 2))+
  geom_shape(data = hull, alpha = 0.7, expand = 0.015, radius = 0.015,
             aes(fill = cluster), color = "white", show.legend = FALSE) +
  geom_point(data = x_filtered, aes(shape = national2), color = "white", size = 2) +
  geom_point(data = x_filtered, aes(color = morphid2, shape = national2)) +
  geom_text_repel(data = labels, aes(label = cluster), size = 4, nudge_y=0.001, nudge_x=0.001)

splitstree_plot

# Save the plot
ggsave("LantCama/outputs/LantCama_splitstree_HDBSCAN_clusters_80miss_maf2.pdf",
       splitstree_plot, width = 20, height = 20, units = "cm", dpi = 600)

ggsave("LantCama/outputs/LantCama_splitstree_HDBSCAN_clusters_80miss_maf2.png",
       splitstree_plot, width = 20, height = 20, units = "cm", dpi = 600)


### UPGMA ###
# write input 
# for iqtree and mrbayes
dart2svdquartets(dms_maf2, RandRbase, species, dataset, add_pop=TRUE, pop=dms_maf2$sample_names)

# read what we just wrote
genotype_matrix <- read.nexus.data('/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen/raw_SNPFilt_1SNPperClone/svdq/LantCama_DLan23-8067.nex')

x <- nexus2DNAbin(genotype_matrix)
x2 <- as.phyDat(x)

dm <- dist.hamming(x2)
tree <- upgma(dm)
# NJ# NJupgma()
set.seed(123)
UPGMAtrees <- bootstrap.phyDat(x2,
                               FUN=function(x)upgma(dist.hamming(x)), bs=1000)
treeUPGMA <- plotBS(tree, UPGMAtrees, "phylogram")

write.tree(treeUPGMA, file='treeUPGMA.tree')

#### read in tree ####

upgma <- read.tree('/Users/eilishmcmaster/Documents/LantCama/upgma/treeUPGMA.tree')
upgma$node.label <- round(as.numeric(upgma$node.label),0)

x1 <- m2[m2$sample %in% upgma$tip.label,]
rownames(x1) <- x1$sample
nation_colours <- named_list_maker(x1$national2, 'Paired',11)

#### plot tree ####
##### full tree ####


ggtree_obj <- ggtree(upgma, size=0.3) %<+% x1 

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('clusters','morphid2','national2')]),#'svdq_pop_label',
                offset=0.0007, width=.05, font.size=0,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("HBDSCAN cluster",'Morphotype',"Country"), hjust=0) +
  scale_fill_manual(values=c(tsne_cols,nation_colours, morphid_colours), na.value = "white") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size=0.8) +
  geom_rootedge(0.0005, size=0.3)+
  # geom_label2(data=upgma, aes(label=node),color='red', nudge_x = 0.0001, label.size=0, fill="transparent", size=1)+
  # geom_label2(data=upgma, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
  #             color='red', nudge_x = 0.0004, label.size=0, fill="transparent", size=2) +
  theme(legend.position = "none",plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)


custom_legend_theme <-   theme(legend.key.size = unit(0.5, 'lines'),
                               legend.title = element_text(size = 8),
                               legend.text = element_text(size = 6))
# Create separate ggplots for each fill scheme
plot_svdq_pop <- ggplot(m2, aes(x=long,y=lat, fill=clusters)) +
  geom_tile() +
  scale_fill_manual(name = "HBDSCAN\ncluster", values = tsne_cols2, na.translate = FALSE) +
  custom_legend_theme

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = nation_colours, na.value = "grey90")+
  custom_legend_theme


plot_morphid2 <- ggplot(m2, aes(x=long,y=lat, fill=morphid2)) +
  geom_tile() +
  scale_fill_manual(name = "Morphotype", values = morphid_colours, na.value = "grey90")+
  custom_legend_theme

plot_prob <- ggplot(data.frame(label=upgma$node.label, x=1:length(upgma$node.label), y=1:length(upgma$node.label)),
                    aes(x=x, y=y,colour=label)) +
  geom_point() +scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)+labs(color="Probability")+custom_legend_theme


# Extract legends from the ggplots
legend_svdq_pop <- cowplot::get_legend(plot_svdq_pop)
legend_national2 <- cowplot::get_legend(plot_national2)
legend_prob <- cowplot::get_legend(plot_prob)
legend_morphid2 <- cowplot::get_legend(plot_morphid2)



legends <- cowplot::plot_grid(legend_prob, legend_svdq_pop,legend_morphid2, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 5/1) # Adjust aspect 
legends2 <- cowplot::plot_grid(legend_morphid2, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 3/1) # Adjust aspect 


combined_plot <- cowplot::plot_grid(hmt, legends,nrow = 1, rel_widths = c(1, 0.15))
combined_plot

ggsave("LantCama/outputs/LantCama_upgma_maf2.pdf",
       combined_plot, width = 20, height = 40, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_upgma_maf2.png",
       combined_plot, width = 25, height = 40, units = "cm", dpi=300)
