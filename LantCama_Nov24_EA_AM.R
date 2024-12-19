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

table(m2$national2)

#### read clusters ###

m_clusters <- read.xlsx('LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx')
m2$clusters <- m_clusters$cluster[match(m2$sample, m_clusters$sample)] %>% as.vector()

#### colours 

countries <- m2 %>%
  filter(national2 == "Native range") %>%
  group_by(country) %>%
  summarise(avg_lat = mean(lat, na.rm = TRUE)) %>%
  arrange(avg_lat) %>%
  pull(country)

m2$country <- factor(m2$country, levels=countries)

country_colours <- named_list_maker(countries, 'RdYlGn', 11)

morphid_colours <- c(pink="#EE6677", PER="forestgreen", red="red3", white="#66CCEE", orange="orange", undetermined="#2B2B2B")

tsne_cols <- structure(c("white", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                         "#FB9A99", "#E31A1C", "#FDBF6F"), names = c(NA, "A", "B", "C", 
                                                                     "D", "E", "F", "G"))
nation_colours <- named_list_maker(m2$national2, 'Paired',11)


tsne_cols2 <- tsne_cols[!is.na(names(tsne_cols))]


custom_theme <- theme(axis.text = element_text(size=8),
                      axis.title = element_text(size=10),
                      legend.text = element_text(size=8),
                      legend.title = element_text(size=10),
                      plot.title = element_text(size = 10),
                      legend.key.size = unit(0.5, 'lines'),
                      legend.key.height = unit(0.5, 'lines'))


custom_legend_theme <-   theme(legend.key.size = unit(0.5, 'lines'),
                               legend.title = element_text(size = 10),
                               legend.text = element_text(size = 8))

#### splitstree ####
# dms_maf2 <- remove.by.maf(dms, 0.02)
# length(dms_maf2$locus_names)
# splitstree(dist(dms_maf2$gt, method = "euclidean"), 'LantCama/outputs/20241205_all_LantCama_nexus_file_for_R_80missing_maf2.nex')

library(tanggle)
library(RSplitsTree)
library(ggplot2)
library(phangorn)
library(ggforce)

# Read network data from Nexus file
Nnet <- phangorn::read.nexus.networx('LantCama/outputs/20241205_all_LantCama_nexus_file_for_R_80missing_maf2.nex')


splitstree <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

splitstree[Nnet$translate$node,"sample"] <- Nnet$translate$label
splitstree <- merge(splitstree, m2, by="sample", all.x=TRUE, all.y=FALSE)
splitstree$cluster[!is.na(splitstree$sample)&is.na(splitstree$cluster)] <- "ungrouped"
# 
net_x_axis <- max(splitstree$x)-min(splitstree$x)
net_y_axis <- max(splitstree$y)-min(splitstree$y)


levels2 <- unique(dms$meta$analyses[,"national2"]) %>% sort()
levels2_shapes <- setNames(1:nlevels(factor(levels2)), levels2)

hull <- splitstree %>% group_by(cluster) %>%
  slice(chull(x, y))
hull <- hull[!hull$cluster=='ungrouped',]


#


# Create hull data
hull <- splitstree %>%
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
  # scale_fill_manual(values = tsne_cols, na.translate = FALSE) +
  scale_colour_manual(values = morphid_colours, na.translate = FALSE) +
  scale_shape_manual(values=16:17,na.translate = FALSE)+
  theme_void() +
  # expand_limits(x = c(min(splitstree$x) - 0.01 * net_x_axis, max(splitstree$x) + 0.1 * net_x_axis),
  #               y = c(min(splitstree$y) - 0.01 * net_y_axis, max(splitstree$y) + 0.01 * net_y_axis)) +
  theme(legend.position = "right", legend.key = element_blank(),legend.justification = "top",
        legend.key.size = unit(0, 'lines')) +#, plot.margin = margin(5,5,5,5)
  coord_fixed() +
  labs(color = "Morphotype", shape = "Origin", fill = "HDBSCAN clusters") +
  guides(colour = guide_legend(title.position = "top", ncol = 1, override.aes = list(fill = NA, linetype = 0)), 
         # fill = guide_legend(title.position = "top", nrow = 2, override.aes = list(color = NA)),
         shape = guide_legend(title.position = "top", nrow = 2))+
  geom_shape(data = hull, alpha = 0.3, expand = 0.015, radius = 0.015,
             aes(group=cluster),fill='black', color = "white", show.legend = FALSE) + #aes(fill = cluster)
  geom_point(data = splitstree, aes(shape = national2), color = "white", size = 2) +
  geom_point(data = splitstree, aes(color = morphid2, shape = national2)) +
  geom_text_repel(data = labels, aes(label = cluster), size = 4, nudge_y=0.001, nudge_x=0.001)+
  custom_legend_theme

splitstree_plot

# Save the plot
ggsave("LantCama/outputs/LantCama_splitstree_HDBSCAN_clusters_80miss_maf2.pdf",
       splitstree_plot, width = 20, height = 20, units = "cm", dpi = 600)

ggsave("LantCama/outputs/LantCama_splitstree_HDBSCAN_clusters_80miss_maf2.png",
       splitstree_plot, width = 20, height = 20, units = "cm", dpi = 600)

# 
# ### UPGMA #################
# # write input
# # for iqtree and mrbayes
# dart2svdquartets(dms_maf2, RandRbase, species, dataset, add_pop=TRUE, pop=dms_maf2$sample_names)
# 
# # read what we just wrote
# genotype_matrix <- read.nexus.data('/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen/raw_SNPFilt_1SNPperClone/svdq/LantCama_DLan23-8067.nex')
# 
# x <- nexus2DNAbin(genotype_matrix)
# x2 <- as.phyDat(x)
# 
# dm <- dist.hamming(x2)
# tree <- upgma(dm)
# # NJ# NJupgma()
# set.seed(123)
# UPGMAtrees <- bootstrap.phyDat(x2,
#                                FUN=function(x)upgma(dist.hamming(x)), bs=10000)
# treeUPGMA <- plotBS(tree, UPGMAtrees, "phylogram")
# 
# write.tree(treeUPGMA, file='treeUPGMA.tree')

#### legends ####

# Create separate ggplots for each fill scheme
plot_cluster_pop <- ggplot(m2, aes(x=long,y=lat, fill=clusters)) +
  geom_tile() +
  scale_fill_manual(name = "HBDSCAN\ncluster", values = tsne_cols2, na.translate = FALSE) +
  custom_legend_theme

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Origin", values = nation_colours, na.value = "grey90")+
  custom_legend_theme


plot_morphid2 <- ggplot(m2, aes(x=long,y=lat, fill=morphid2)) +
  geom_tile() +
  scale_fill_manual(name = "Morphotype", values = morphid_colours, na.value = "grey90")+
  custom_legend_theme

plot_country <- ggplot(m2[which(m2$national2=="Native range"& !is.na(m2$country)),], aes(x=long,y=lat, fill=country)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = country_colours)+
  custom_legend_theme

plot_prob <- ggplot(data.frame(label=upgma$node.label, x=1:length(upgma$node.label), y=1:length(upgma$node.label)),
                    aes(x=x, y=y,colour=label)) +
  geom_point() +scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)+labs(color="Probability")+custom_legend_theme


# Extract legends from the ggplots
legend_cluster_pop <- cowplot::get_legend(plot_cluster_pop +theme(legend.direction = 'horizontal'))
legend_national2 <- cowplot::get_legend(plot_national2+theme(legend.direction = 'horizontal'))
legend_prob <- cowplot::get_legend(plot_prob)
legend_morphid2 <- cowplot::get_legend(plot_morphid2+theme(legend.direction = 'horizontal'))
legend_country <- cowplot::get_legend(plot_country+theme(legend.direction = 'horizontal'))


legends <- cowplot::plot_grid(legend_cluster_pop,legend_morphid2, legend_national2, nrow=1, align="hv")  # Adjust aspect #legend_prob
legends2 <- cowplot::plot_grid(legend_morphid2, legend_national2,  nrow=1, align="hv") #+ theme(aspect.ratio = 1/3) # Adjust aspect 
legends3 <- cowplot::plot_grid(legends2, legend_country,nrow=2)

#### read in tree ####

upgma <- read.tree('/Users/eilishmcmaster/Documents/LantCama/upgma/treeUPGMA.tree')

upgma$node.label <- round(as.numeric(upgma$node.label),0)
# # upgma <- unroot(upgma)
# upgma <- ape::root(upgma, node = 990, resolve.root=TRUE)

x1 <- m2[m2$sample %in% upgma$tip.label,]
rownames(x1) <- x1$sample

#### plot tree ####
##### full tree ####

ggtree_obj <- ggtree(upgma, size=0.3) %<+% x1 

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('clusters','morphid2','national2', 'country')]),#'svdq_pop_label',
                offset=0.0007, width=.05, font.size=0,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("HBDSCAN cluster",'Morphotype',"Origin",'Country'), hjust=0) +
  scale_fill_manual(values=c(tsne_cols,nation_colours, morphid_colours, country_colours), na.value = "white") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size=0.8, align = TRUE, linetype = "dotted", linesize=0.1) +
  geom_rootedge(0.0005, size=0.3)+
  geom_label2(data=upgma, aes(label=node),color='red', nudge_x = 0.0001, label.size=0, fill="transparent", size=1)+
  # geom_label2(data=upgma, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
  #             color='red', nudge_x = 0.0003, label.size=0, fill="transparent", size=2) +
  theme(legend.position = "none",plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))
  # geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)


combined_plot <- cowplot::plot_grid(hmt, legends3 ,nrow = 2, rel_heights = c(1, 0.15))

ggsave("LantCama/outputs/LantCama_upgma_maf2.pdf",
       combined_plot, width = 25, height = 40, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_upgma_maf2.png",
       combined_plot, width = 25, height = 40, units = "cm", dpi=300)

#### collapsed UPGMA ####
library(ggtree)
library(ggplot2)
library(ape)
library(openxlsx)
library(tidytree)
library(RColorBrewer)

nodes.to.collapse <- c(1293,1281,1003, 1045, 946, 665, 784)
nodes.identity <- c('G','E', 'F','B', 'C', 'A', 'D')

names(nodes.identity) <- nodes.to.collapse

cols_for_clades <- c(A="#EE6677", B="forestgreen", C="red3", D="#66CCEE", E="orange", `F`="#EE6677", G='orange')

# Scale clade function proportional to the number of samples
scaleMyClade <- function(.p, .node) {
  offs <- offspring(.p$data, .node)
  num_samples <- nrow(offs)
  scaling_factor <- 1 / (log10(num_samples) + 1) # Adjust the scaling factor as needed
  scaleClade(.p, .node, scaling_factor)
}


# Collapse clade function with labels and colors
collapseMyClade <- function(.p, .node) {
  label <- nodes.identity[as.character(.node)]
  fill_color <- cols_for_clades[label]
  .p <- collapse(.p, .node, "max",  fill = fill_color, color="black", size=0.3)
  .p <- .p + geom_cladelabel(node = .node,align = TRUE, vjust=-0.02,offset.text=-0.0005, label = label,
                             fontsize = 3)
  return(.p)
}


# Apply scaling and collapsing to the nodes
ggtree_obj2 <- Reduce(scaleMyClade, nodes.to.collapse, ggtree_obj)
ggtree_obj2 <- Reduce(collapseMyClade, nodes.to.collapse, ggtree_obj2)

gheatmap(ggtree_obj2, as.matrix(x1[,c('morphid2','national2')]),
         offset = 0.001, width = .1, font.size = 2,
         colnames_position = "top")
         # custom_column_labels = c('Morphotype', "Origin"))
  # theme(legend.position = "none", plot.margin = margin(3, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  # coord_cartesian(clip = "off")


# Add gheatmap and additional layers as required
hmt2 <- gheatmap(ggtree_obj2, as.matrix(x1[,c('morphid2','national2','country')]),#'svdq_pop_label',
                 offset = 0.001, width = .1, font.size = 10,
                 colnames_angle = 90, colnames_position = "top",
                 custom_column_labels = c('Morphotype', "Origin", 'Country'), hjust = 0) +
  scale_fill_manual(values = c(nation_colours, morphid_colours,country_colours), na.value = "white") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size = 0.8) +
  geom_rootedge(0.0005, size = 0.3) +
  theme(legend.position = "none", plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
              color='red', nudge_x = -0.0003,nudge_y=1.1, label.size=0, fill="transparent", size=1.5)
# geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)


combined_plot2 <- cowplot::plot_grid(hmt2, legends3, nrow = 2, rel_heights = c(1, 0.12))

g1 <- combined_plot2
g2 <- splitstree_plot
# Convert g2 into a grob
g2_grob <- ggplotGrob(g2)

# Overlay g2 on g1 using annotation_custom
g1_with_g2 <- g1 +
  annotation_custom(
    grob = g2_grob,
    xmin = 0, xmax = 0.65, # Adjust coordinates as needed
    ymin = 0.58, ymax = 0.99   # Adjust coordinates as needed
  )

# Print the combined plot
# g1_with_g2

ggsave("LantCama/outputs/LantCama_upgma_splitstree_maf2.pdf",
       g1_with_g2, width = 20, height = 25, units = "cm", dpi=600, bg='white')

ggsave("LantCama/outputs/LantCama_upgma_splitstree_maf2.png",
       g1_with_g2, width = 18, height = 24, units = "cm", dpi=300, bg='white')

######## heterozygosity ##########
meta <- read.xlsx("/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx") # get meta for lat long limits

ho_0 <- individual_Ho(dms$gt, genetic_group_variable = rep("lcamara",length(dms$sample_names)), maf = 0)
ho_2 <- individual_Ho(dms$gt, genetic_group_variable = rep("lcamara",length(dms$sample_names)), maf = 0.02)

hist(ho_2)
m2$Ho_0 <- ho_0[match(m2$sample, names(ho_0))]
m2$Ho_2 <- ho_2[match(m2$sample, names(ho_2))]
m2$cluster <- meta$cluster[match(m2$sample, meta$sample)]
m2$cluster[which(m2$national2=="Native range")] <- "Native range"
m2$cluster[which(is.na(m2$cluster))] <- "Unclustered Australian"
m2$cluster <- factor(m2$cluster, levels = c('A','B','C','D','E','F','G','Unclustered Australian','Native range'))
tsne_cols3 <- c(tsne_cols,`Unclustered Australian`='grey80', `Native range`='grey30')

ho_boxplot <- ggplot(m2, aes(x=cluster, y=Ho_2, fill=cluster))+
  geom_boxplot()+
  scale_fill_manual(values=tsne_cols3)+
  theme_few()+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, hjust=1))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 6))+
  labs(x="HDBSCAN Cluster", y="Ho (MAF 2%)")
ho_boxplot


ho_hist <- ggplot(m2, aes(x=Ho_2, fill=cluster))+
  geom_histogram(position='stack', binwidth = 0.01, boundary=0)+
  scale_fill_manual(values=tsne_cols3)+
  theme_few()+
  theme(legend.position = 'bottom', 
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(0.5, 'lines')
        )+
  labs(x="Ho (MAF 2%)", y="Count (individuals)", fill="HDBSCAN\nCluster")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits = c(0,170))+
  guides(fill=guide_legend(ncol=3))

ho_hist

combined_hos <- ggarrange(ho_boxplot, ho_hist, labels="auto", font.label = list(face = "plain"),align='v')
ggsave('LantCama/outputs/Figure2_ho_box_hist.png',dpi = 600, combined_hos, width = 20, height = 10, units = "cm")


#### range distances ####

library(dplyr)
library(geosphere)  # For geodistances

# Define a function to calculate the max distance in a cluster
max_distance <- function(lat, long) {
  coords <- cbind(long, lat)  # Longitude first for geosphere functions
  coords <- na.omit(coords)  # Remove rows with NA values
  distances <- distm(coords)  # Calculate pairwise distances
  max(distances)  # Return the maximum distance
}

# Group by cluster and calculate the maximum distance
cluster_max_distances <- meta %>%
  group_by(cluster) %>%
  summarise(
    max_dist_km = max_distance(lat, long) / 1000  # Convert meters to kilometers
  )

# View the result
print(cluster_max_distances)

# Calculate the latitudinal span for each cluster
latitudinal_span <- meta %>%
  group_by(cluster) %>%
  summarise(
    min_lat = min(lat, na.rm = TRUE),  # Minimum latitude, ignoring NA
    max_lat = max(lat, na.rm = TRUE),  # Maximum latitude, ignoring NA
    degrees_latitude = max_lat - min_lat  # Difference in degrees
  )

# View the result
print(latitudinal_span)
