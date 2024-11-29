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
d3        <- dart.meta.data.merge(d1, meta) %>% remove.by.list(.,meta$sample_names[which(!is.na(meta$analyses[,'EA_only']))])
d4        <- remove.poor.quality.snps(d3, min_repro=0.98, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
length(d5$locus_names)
dms <- remove.by.missingness(d5, 0.8)
length(dms$sample_names)

m2 <- dms$meta$analyses %>% as.data.frame
m2$lat <- as.numeric(m2$lat)
m2$long <- as.numeric(m2$long)
#### colours 
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")
svdq_pop_colours <- named_list_maker(m2$svdq_pop, 'Spectral', 11)
svdq_pop_colours <- c(svdq_pop_colours, 'ungrouped'='grey30')


### PCA  ################################################################################
dms_maf2 <- remove.by.maf(dms, 0.02)

gen_d5 <- new("genlight", dms_maf2[["gt"]]) #convert df to genlight object for glPca function
gen_pca <- glPca(gen_d5, parallel=TRUE, nf=6) #do pca -- this method somehow allows the input to hav1 NAs

g_pca_df <- gen_pca[["scores"]] #extract PCs
g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # some in DArT are not in meta?

pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:6]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes

pca_plot1 <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=morphid2))+ xlab(pcnames[1])+ylab(pcnames[2])+
  geom_point(size=0.5)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="Morphotype")+
  theme(legend.key.size = unit(0, 'lines'))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)
pca_plot1

pca_plot2 <- ggplot(g_pca_df2, aes(x=PC3, y=PC4, colour=morphid2))+ xlab(pcnames[3])+ylab(pcnames[4])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top")) +
  scale_colour_manual(values=morphid_colours)

pca_plot3 <- ggplot(g_pca_df2, aes(x=PC5, y=PC6, colour=morphid2))+ xlab(pcnames[5])+ylab(pcnames[6])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)


all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3, labels=c("A","B","C"),
                            common.legend = TRUE, ncol=3, legend = "bottom")
all3_pca_plots


### TSNE ################################################################################
# Existing code to compute the distance matrix
genotype_matrix <- dms$gt %>% as.matrix()
d <- dist(genotype_matrix, method = "euclidean")
d_matrix <- as.matrix(d)

# Perform t-SNE on the distance matrix (or you could use the original genotype matrix)
tsne_result <- Rtsne(d, dims = 2, perplexity = 15, 
                     check_duplicates = FALSE, theta=0,
                     is_distance=TRUE,
                     num_threads=4)

# Run HDBSCAN on the t-SNE output (2D coordinates)
tsne_df <- data.frame(tsne_result$Y)  # Extract the 2D coordinates from t-SNE
rownames(tsne_df) <- names(d)
colnames(tsne_df) <- c("tSNE1", "tSNE2")  # Rename columns for clarity

# Apply HDBSCAN to the t-SNE coordinates (2D space)
hdbscan_result <- hdbscan(tsne_df, minPts = 5)

hdbscan_clusters <- hdbscan_result$cluster

#### determine monophyletic clusters ####
# Create a function to check if a given cluster is monophyletic
is_monophyletic_cluster <- function(tree, cluster_assignment, cluster_id) {
  # Find the indices of tips that belong to the cluster
  cluster_tips <- which(cluster_assignment == cluster_id)
  # Check if the cluster_tips form a monophyletic clade
  return(is.monophyletic(tree, cluster_tips))
}

# Convert hclust object to phylo tree
phylo_tree <- as.phylo(hclust(d, method='average'))

# Check monophyly for each unique cluster
cluster_ids <- unique(hdbscan_clusters)
monophyly_results <- sapply(cluster_ids, function(cluster_id) {
  is_monophyletic_cluster(phylo_tree, hdbscan_clusters, cluster_id)
})

monophyletic_clusters <- cluster_ids[monophyly_results]

#### filter clusters ####

min_cluster_size <- 10
hdb_df2 <- data.frame(sample=names(d),hdb_cluster=hdbscan_result$cluster)

small_clusters <- names(which(table(hdb_df2$hdb_cluster) < min_cluster_size))
# Set the cluster value to 0 for those clusters
hdb_df2$hdb_cluster[hdb_df2$hdb_cluster %in% small_clusters] <- 0
# hdb_df2$hdb_cluster[hdb_df2$sample %in% samples_to_remove2] <- 0
hdb_df2$hdb_cluster[!(hdb_df2$hdb_cluster %in% monophyletic_clusters)] <- 0

hdb_df2$hdb_cluster[which(hdb_df2$hdb_cluster==0)] <- NA

hdb_df2 <- merge(hdb_df2, m2, by='sample')
hdb_df2 <- merge(hdb_df2, tsne_df,by.x='sample', by.y=0)
hdb_df2 <- hdb_df2 %>% arrange(lat)
hdb_cluster_unique <- unique(hdb_df2$hdb_cluster)
hdb_cluster_unique <- hdb_cluster_unique[!is.na(hdb_cluster_unique)]
new_cluster_alphabetical <- LETTERS[1:length(hdb_cluster_unique)]
hdb_df2$cluster <- new_cluster_alphabetical[match(hdb_df2$hdb_cluster, hdb_cluster_unique)]

tsne_cols <- c(rainbow(length(unique(hdb_df2$cluster))))
names(tsne_cols) <- unique(hdb_df2$cluster)
# Plot the t-SNE result with HDBSCAN clusters

hull <- hdb_df2 %>% 
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>% 
  slice(chull(tSNE1, tSNE2))

tsne_plot1 <- ggplot(hdb_df2, aes(x = tSNE1, y = tSNE2, color = factor(cluster))) +
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  geom_shape(data=hull, mapping=aes(group=cluster), color="black", linetype='dashed', alpha=0.1,
             expand = 0.02, radius=0.01)+
  geom_point(size=0.5) +
  theme(legend.key.size = unit(0, 'lines'))+
  scale_color_manual(values = tsne_cols) +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "HDBSCAN Cluster") 

tsne_plot2 <- ggplot(hdb_df2, aes(x = tSNE1, y = tSNE2, color = morphid2)) +
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  geom_shape(data=hull, mapping=aes(group=cluster), color="black", linetype='dashed', alpha=0.1,
             expand = 0.02, radius=0.01)+
  theme(legend.key.size = unit(0, 'lines'))+
  geom_point(size=0.5) +
  scale_color_manual(values = morphid_colours) +
  labs( x = "t-SNE 1", y = "t-SNE 2", color = "Morphotype")

tsne_plot1

ggarrange(pca_plot1, tsne_plot2, tsne_plot1, nrow=3, align="hv")

#### dist heatmap ####
d_matrix2 <- merge(d_matrix, hdb_df2, by.x=0, by.y='sample')
rownames(d_matrix2) <- d_matrix2$Row.names
d_matrix2$Row.names <- NULL
d_matrix2 <- d_matrix2[match(colnames(d_matrix2)[1:nrow(d_matrix2)],rownames(d_matrix2)),]

tsne_cols2 <- tsne_cols[!is.na(names(tsne_cols))]

# Create row and column annotations
row_anno <- rowAnnotation(
  cluster = as.factor(d_matrix2$cluster),
  col = list(cluster = tsne_cols2),
  annotation_legend_param = list(
    cluster = list(title = "HDBSCAN Cluster")
  ),
  annotation_name_gp = gpar(fontsize = 0),
  na_col = "white"
  )

row_anno2 <- rowAnnotation(
  Morphotype = d_matrix2$morphid2,
  col = list(Morphotype = morphid_colours),
  annotation_legend_param = list(
    Morphotype = list(title = "Morphotype")
  ),
  annotation_name_gp = gpar(fontsize = 0)
)

row_anno3 <- rowAnnotation(
  Latitude = as.numeric(d_matrix2$lat),
  col = list(
    Latitude = colorRamp2(
      c(min(d_matrix2$lat, na.rm = TRUE), mean(d_matrix2$lat, na.rm = TRUE), max(d_matrix2$lat, na.rm = TRUE)),
      c("red", "white", "blue")
    )
  ),
  annotation_name_gp = gpar(fontsize = 0)
)

col_anno <- HeatmapAnnotation(
  cluster = as.factor(d_matrix2$cluster),
  col = list(cluster = tsne_cols2),
  annotation_legend_param = list(
    cluster = list(title = "HDBSCAN Cluster")
  ),
  annotation_name_gp = gpar(fontsize = 0),
  na_col = "white"
  )

# Create the heatmap
ht <- Heatmap(
  d_matrix2[,1:nrow(d_matrix2)] %>% as.matrix,
  name = "Distance",
  clustering_method_columns = "average",
  clustering_method_rows = "average",
  clustering_distance_rows=d,
  clustering_distance_columns = d,
  top_annotation = col_anno,
  right_annotation = c(row_anno, row_anno2,row_anno3),
  show_column_names = FALSE,  # Remove column names
  show_row_names = FALSE,   
  col = colorRamp2(c(min(d_matrix), max(d_matrix)), c("white", "black")),
  row_dend_width = unit(3, "cm"),          # Adjust row dendrogram size
  column_dend_height = unit(3, "cm")
)

draw(ht, merge_legends = TRUE)

combined_plots <- multi_panel_figure(
  width = c(7.5, 25),   # Adjust these dimensions as needed
  height = c(7,7,7),
  unit = "cm",
  panel_label_type = "upper-roman"
)

# Fill the panels with the respective plots
combined_plots %<>%
  fill_panel(pca_plot1 + theme(legend.position = "none"), column = 1, row = 1, label = "A") %<>%
  fill_panel(tsne_plot2+ theme(legend.position = "none"), column = 1, row = 2, label = "B") %<>%
  fill_panel(tsne_plot1+ theme(legend.position = "none"), column = 1, row = 3, label = "C") %<>%
  fill_panel(draw(ht, merge_legends = TRUE), column = 2, row = 1:3, label = "D")

ggsave('LantCama/outputs/Figure1_combined_plots.pdf', combined_plots, width = 34, height = 23, units = "cm")


# ##### LEA ####
# 
# library(LEA)
# 
# nd_lea <- dart2lea(dms, RandRbase, species, dataset)
# kvalrange <- 1:15
# # 
# # snmf1 <- snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 10, project = "new", CPU=8)
# # 
# # save(snmf1, file='LantCama/popgen/LantCama_EA_only_snmf.RData')
# 
load(file='LantCama/popgen/LantCama_EA_only_snmf.RData')

plot(snmf1, col = "blue", pch = 19, cex = 1.2)
best = which.min(cross.entropy(snmf1, K = K_chosen))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold",'blue')

K_chosen <- 12
# 
# barchart(snmf1, K = K_chosen, run = best,
#          border = NA, space = 0,
#          col = my.colors,
#          xlab = "Individuals",
#          ylab = "Ancestry proportions",
#          main = "Ancestry matrix") -> bp
# 
# axis(1, at = 1:length(bp$order),
#      labels = bp$order, las=1,
#      cex.axis = .4)
# 
# qmatrix_df <- as_tibble(Q(snmf1, K = K_chosen, run=which.min(cross.entropy(snmf1, K = K_chosen)))) %>%
#   mutate(sample = dms$sample_names) %>%
#   pivot_longer(-sample, names_to = "lea_cluster", values_to = "proportion")
# 
# qmatrix_df2 <- merge(qmatrix_df, hdb_df2, by='sample')
# 
# ggplot(qmatrix_df2, aes(x = sample, y = proportion, fill = lea_cluster)) +
#   geom_bar(stat = "identity", width = 1) +
#   facet_grid(~cluster, scales = "free_x", space="free_x")+
#   theme_few() +
#   scale_fill_brewer(palette = "Set3") + # Use a more subtle color palette
#   theme(axis.text.x = element_blank())+
#   scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
#   labs(x = "Individuals", y = "Ancestry Proportion", fill = "Cluster")
