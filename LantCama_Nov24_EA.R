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

EA_only <- c(meta$sample_names[which(!is.na(meta$analyses[,'EA_only']))])

d3.1 <- remove.by.list(d3,EA_only)
d4        <- remove.poor.quality.snps(d3.1, min_repro=0.98, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
length(d5$locus_names)
dms <- remove.by.missingness(d5, 0.8)
length(dms$sample_names)

m2 <- dms$meta$analyses %>% as.data.frame
m2$lat <- as.numeric(m2$lat)
m2$long <- as.numeric(m2$long)
#### colours 
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")


### PCA  ################################################################################
dms_maf2 <- remove.by.maf(dms, 0.02)
length(dms_maf2$locus_names)

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
  labs(colour="`HBDSCAN Cluster`")+
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
tsne_result <- Rtsne(d, dims = 2, perplexity = 18, 
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

min_cluster_size <- 11
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



n_clusters <- length(unique(hdb_df2$cluster[!is.na(hdb_df2$cluster)]))  # Exclude NA from the count
# tsne_cols <- brewer.pal(n_clusters, "Paired")
# tsne_cols <- c("white", tsne_cols)
# names(tsne_cols) <- c(NA, unique(hdb_df2$cluster[!is.na(hdb_df2$cluster)]))

tsne_cols <- structure(c("white", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                         "#FB9A99", "#E31A1C", "#FDBF6F"), names = c(NA, "A", "B", "C", 
                                                                     "D", "E", "F", "G"))

#### write_out_clusters ####
out_cols <- c('sample','lat','long','site','site_original','morphid2','cluster')
write.xlsx(hdb_df2[,out_cols], 'LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx')


#### plot tsne ####
# Plot the t-SNE result with HDBSCAN clusters
hull <- hdb_df2 %>% 
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>% 
  slice(chull(tSNE1, tSNE2))

tsne_plot1 <- ggplot(hdb_df2, aes(x = tSNE1, y = tSNE2, color = factor(cluster))) +
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  geom_shape(data=hull, mapping=aes(group=cluster), color="black", linetype='dashed', alpha=0,
             expand = 0.02, radius=0.01)+
  geom_point(size=0.5) +
  theme(legend.key.size = unit(0, 'lines'))+
  scale_color_manual(values = tsne_cols) +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "HDBSCAN Cluster") 

tsne_plot2 <- ggplot(hdb_df2, aes(x = tSNE1, y = tSNE2, color = morphid2)) +
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  geom_shape(data=hull, mapping=aes(group=cluster), color="black", linetype='dashed', alpha=0,
             expand = 0.02, radius=0.01)+
  theme(legend.key.size = unit(0, 'lines'))+
  geom_point(size=0.5) +
  scale_color_manual(values = morphid_colours) +
  labs( x = "t-SNE 1", y = "t-SNE 2", color = "`HBDSCAN Cluster`")

tsne_plot1


#### dist heatmap ####
d_matrix2 <- merge(d_matrix, hdb_df2, by.x=0, by.y='sample')
rownames(d_matrix2) <- d_matrix2$Row.names
d_matrix2$Row.names <- NULL
d_matrix2 <- d_matrix2[match(colnames(d_matrix2)[1:nrow(d_matrix2)],rownames(d_matrix2)),]

tsne_cols2 <- tsne_cols[!is.na(names(tsne_cols))]

# Create row and column annotations
row_anno <- rowAnnotation(
  `HDBSCAN Cluster` = as.factor(d_matrix2$cluster),
  col = list(`HDBSCAN Cluster` = tsne_cols2),
  annotation_legend_param = list(
    `HDBSCAN Cluster` = list(title = "HDBSCAN\nCluster")
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),
  na_col = "white"
  )

row_anno2 <- rowAnnotation(
  Morphotype = d_matrix2$morphid2,
  col = list(Morphotype = morphid_colours),
  annotation_legend_param = list(
    Morphotype = list(title = "Morphotype")
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10)
)

row_anno3 <- rowAnnotation(
  Latitude = as.numeric(d_matrix2$lat),
  col = list(
    Latitude = colorRamp2(
      c(min(d_matrix2$lat, na.rm = TRUE), mean(d_matrix2$lat, na.rm = TRUE), max(d_matrix2$lat, na.rm = TRUE)),
      c("red", "white", "blue")
    )
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10)
)

col_anno <- HeatmapAnnotation(
  cluster = as.factor(d_matrix2$cluster),
  col = list(cluster = tsne_cols2),
  annotation_legend_param = list(
    cluster = list(title = "HDBSCAN\nCluster")
  ),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 0),
  na_col = "white"
  )

# Create the heatmap
ht <- Heatmap(
  d_matrix2[,1:nrow(d_matrix2)] %>% as.matrix,
  name = "Distance",
  show_row_dend = FALSE,
  # show_column_dend = FALSE,
  clustering_method_columns = "average",
  clustering_method_rows = "average",
  clustering_distance_rows=d,
  clustering_distance_columns = d,
  top_annotation = col_anno,
  right_annotation = c(row_anno, row_anno2,row_anno3),
  show_column_names = FALSE,  # Remove column names
  show_row_names = FALSE,   
  show_heatmap_legend = FALSE,
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

# add clusters to DMS
dms$meta$cluster <- hdb_df2$cluster[match(dms$sample_names, hdb_df2$sample)] %>% as.vector()
dms$meta$site_cluster <- paste0(dms$meta$site, ifelse(is.na(dms$meta$cluster), "", paste0("(",dms$meta$cluster,")")))

m2$cluster <- hdb_df2$cluster[match(m2$sample,hdb_df2$sample)] %>% as.vector()
m2$site_cluster <- paste0(m2$site, ifelse(is.na(m2$cluster), "", paste0("(",m2$cluster,")")))
# ### LEA ####
# #
library(LEA)
# 
# nd_lea <- dart2lea(dms, RandRbase, species, dataset)
# kvalrange <- 1:20
# snmf1 <- snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 3, project = "new", CPU=8)
# 
# save(snmf1, file='LantCama/popgen/LantCama_EA_only_snmf.RData')
# # #
load(file='LantCama/popgen/LantCama_EA_only_snmf.RData')

K_chosen <- 5
best = which.min(cross.entropy(snmf1, K = K_chosen))
plot(snmf1, col = "blue", pch = 19, cex = 1.2)


qmatrix_df <- as_tibble(Q(snmf1, K = K_chosen, run=which.min(cross.entropy(snmf1, K = K_chosen)))) %>%
  mutate(sample = dms$sample_names) %>%
  pivot_longer(-sample, names_to = "lea_cluster", values_to = "proportion")

qmatrix_df2 <- merge(qmatrix_df, hdb_df2, by='sample')

qmatrix_df2$cluster[is.na(qmatrix_df2$cluster)] <- "Unclustered"

qmatrix_df2 <- qmatrix_df2 %>%
  mutate(lea_cluster = gsub("V", "", lea_cluster))

# Order samples by latitude
qmatrix_df2 <- qmatrix_df2 %>%
  arrange(lat) %>%  # Arrange by latitude
  mutate(sample = factor(sample, levels = unique(sample)))  # Preserve order in factor


#


# Plot ancestry proportions with latitude annotation
main_plot <- ggplot(qmatrix_df2, aes(x = sample, y = proportion, fill = factor(lea_cluster))) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  theme_few() +
  # scale_fill_brewer(palette = "Greys") +  # Use a more subtle color palette
  scale_fill_brewer(palette = "Set3") +  # Use a more subtle color palette
  theme(
    axis.title.y = element_text(size = 11, angle = 0, hjust = 1, vjust = 0.5),  # Right justified and centered vertically
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),  # Remove axis ticks
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(2, 2, 0, 2, unit = "pt")  # Smaller margins (top, right, bottom, left)
  ) +
  scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
  labs(y = "Ancestry Proportion", fill = "")



lat_plot <- ggplot(qmatrix_df2, aes(x = sample, y = 1, fill = lat)) +
  geom_tile() +
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "red", 
    mid = "white", 
    high = "blue", 
    midpoint = mean(qmatrix_df2$lat, na.rm = TRUE),
    name = "Latitude"
  ) +
  labs(y = "Latitude") +
  theme_void() +
  scale_y_continuous(limits = c(0.5, 1.501), expand = c(0, 0)) +
  theme(
    axis.title.y = element_text(size = 11, angle = 0, hjust = 1, vjust = 0.5),  # Add this to show the y-axis label
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # Add panel border
    strip.text = element_blank()  # Remove facet panel text
    # plot.margin = margin(0, 2, 2, 2, unit = "pt")  # Smaller margins (top, right, bottom, left)
  )


morphotype_plot <- ggplot(qmatrix_df2, aes(x = sample, y = 1, fill = morphid2)) +
  geom_tile() +
  scale_fill_manual(values=morphid_colours, na.value = "white")+
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  labs(fill = "", y="Morphotype") +
  theme_void() +
  scale_y_continuous(limits = c(0.5, 1.501), expand = c(0, 0)) +
  theme(
    axis.title.y = element_text(size = 11, angle = 0, hjust = 1, vjust = 0.5),  # Add this to show the y-axis label
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # Add panel border
    strip.text = element_blank()  # Remove facet panel text
    # plot.margin = margin(0, 2, 0, 2, unit = "pt")  # Smaller margins (top, right, bottom, left)
  )

cluster_plot <- ggplot(qmatrix_df2, aes(x = sample, y = 1, fill = cluster)) +
  geom_tile() +
  scale_fill_manual(values=tsne_cols, na.value = "white")+
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  labs(fill = "Latitude", y="HDBSCAN Cluster") +
  theme_void() +
  scale_y_continuous(limits = c(0.5, 1.501), expand = c(0, 0)) +
  theme(
    axis.title.y = element_text(size = 11, angle = 0, hjust = 1, vjust = 0.5),  # Add this to show the y-axis label
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # Add panel border
    strip.text = element_blank()  # Remove facet panel text
    # plot.margin = margin(0, 2, 0, 2, unit = "pt")  # Smaller margins (top, right, bottom, left)
  )



combined_lea_plot <- ggarrange(main_plot,cluster_plot,morphotype_plot,lat_plot,
          nrow=4,
          heights=c(1,0.1,0.1, 0.1),
          align="v", common.legend = TRUE, legend='right')

# combined_lea_plot <-ggarrange(main_plot,morphotype_plot,lat_plot,
#                               nrow=4, 
#                               heights=c(1,0.15,0.15, 0.15),
#                               align="v", common.legend = TRUE, legend='right')

# combined_plots2 <- multi_panel_figure(
#   width = c(7.5, 25),   # Adjust these dimensions as needed
#   height = c(7,7,7,6),
#   unit = "cm",
#   panel_label_type = "upper-roman"
# )
# 
# # Fill the panels with the respective plots
# combined_plots2 %<>%
#   fill_panel(pca_plot1 + theme(legend.position = "none"), column = 1, row = 1, label = "A") %<>%
#   fill_panel(tsne_plot2 + theme(legend.position = "none"), column = 1, row = 2, label = "B") %<>%
#   fill_panel(tsne_plot1 + theme(legend.position = "none"), column = 1, row = 3, label = "C") %<>%
#   fill_panel(draw(ht, merge_legends = TRUE), column = 2, row = 1:3, label = "D") %<>%
#   fill_panel(combined_lea_plot, column = 1:2, row = 4, label = "E")
# 
# 
# ggsave('LantCama/outputs/Figure1_combined_plots2.pdf', combined_plots2, width = 34, height = 30, units = "cm")
# ggsave('LantCama/outputs/Figure1_combined_plots2.png', combined_plots2, width = 34, height = 30,dpi=600, units = "cm")
# 
# ###

combined_plots2 <- multi_panel_figure(
  width = c(6, 18),   # Adjust these dimensions as needed
  height = c(5,5,5,5),
  unit = "cm",
  panel_label_type = "upper-roman"
)

# Fill the panels with the respective plots
combined_plots2 %<>%
  fill_panel(pca_plot1 + theme(legend.position = "none"), column = 1, row = 1, label = "a") %<>%
  fill_panel(tsne_plot2 + theme(legend.position = "none"), column = 1, row = 2, label = "b") %<>%
  fill_panel(tsne_plot1 + theme(legend.position = "none"), column = 1, row = 3, label = "c") %<>%
  fill_panel(draw(ht, merge_legends = TRUE), column = 2, row = 1:3, label = "d") %<>%
  fill_panel(combined_lea_plot, column = 1:2, row = 4, label = "e")


ggsave('LantCama/outputs/Figure1_combined_plots2.pdf', combined_plots2, width = 25, height = 22.1, units = "cm")
ggsave('LantCama/outputs/Figure1_combined_plots2.png', combined_plots2, width = 25, height = 22.1,dpi=600, units = "cm")


### FST ###
# remove site_clusters where n=4
sppop_freq <- as.data.frame(table(dms$meta$site_cluster))
not_n1_site_clusters <- as.vector(sppop_freq[sppop_freq$Freq<5,1]) #remove groups where n<=1
not_n1_samples <- dms$sample_names[which(!(dms$meta$site_cluster %in% not_n1_site_clusters)& !is.na(dms$meta$site_cluster))]
fst_dms <- remove.by.list(dms, not_n1_samples)

length(fst_dms$sample_names)
length(unique(fst_dms$meta$site_cluster))
length(fst_dms$locus_names)

gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$site_cluster, RandRbase,species,dataset, maf_val=0.02, miss_val=0.8) #calculates genetic distance
pS        <- population.pw.spatial.dist(fst_dms, fst_dms$meta$site_cluster) #calculates geographic distance between populations

####plot IBD plot

library(reshape2) #for melting data
library(vegan) #for mantel test

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

#Mantel test
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE) #mantel test, finds if matrices are signficantly similar
man

# mantel plot
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000

# adding metadata for site_clusters
Fst_sig2 <- merge(Fst_sig, distinct(m2[,c("site_cluster","cluster")]), by.x="Var1", by.y="site_cluster", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(m2[,c("site_cluster","cluster")]), by.x="Var2", by.y="site_cluster", all.y=FALSE)
Fst_sig2$same_cluster <- ifelse(Fst_sig2$cluster.x == Fst_sig2$cluster.y, "Intra-cluster", "Inter-cluster")

library(ggforce)
fstp1 <- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_cluster))+
  geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1


ggsave('LantCama/outputs/Figure2_fst.png',dpi = 300, fstp1, width = 12, height = 8, units = "cm")
ggsave('LantCama/outputs/Figure2_fst.pdf',dpi = 300, fstp1, width = 12, height = 8, units = "cm")

paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif)

# Make heatmaps
# geo dist
geo_d <-pS$S #this is a square matrix
mat <- geo_d/1000 # convert to km

#FST
mat2 <-pFst$Fst
diag(mat2) <- NA

order_hm <- Heatmap(mat2,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
od <- colnames(mat2)[column_order(order_hm)]

mat = mat[od, od]
mat2 = mat2[od, od]

# order_hm <- Heatmap(mat,
#                     cluster_rows = TRUE,
#                     cluster_columns = TRUE)
# od <- colnames(mat)[column_order(order_hm)]
#
# mat = mat[od, od]
# mat2 = mat2[od, od]

# agg <- unique(m2[, c("site_cluster", "cluster",'national2')]) # create aggregated df of pop_largeecies and site_cluster
agg <- unique(m2[, c("site_cluster",'cluster')]) # create aggregated df of pop_largeecies and site_cluster

mat2 <- merge(mat2, agg, by.x=0, by.y="site_cluster", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

row_group_ann <- rowAnnotation(`HBDSCAN Cluster` = mat2$cluster,
                               col=list(`HBDSCAN Cluster`=tsne_cols2),
                               na_col="white",
                               annotation_legend_param = list(labels_gp=gpar(fontsize=8),
                                                              title_gp=gpar(fontsize=10)),
                               annotation_name_gp = gpar(fontsize = 0),
                               annotation_name_side="top")



bottom_group_ann <- HeatmapAnnotation(`HBDSCAN Cluster` = mat2$cluster, col = list(`HBDSCAN Cluster` = tsne_cols2),
                                      annotation_name_gp = gpar(fontsize = 0),
                                      annotation_legend_param = list(labels_gp=gpar(fontsize=8),
                                                                     title_gp=gpar(fontsize=10)),
                                      annotation_name_side="left",
                                      na_col = "white")

# specify fst heatmap colours
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))


#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat, na.rm=TRUE)), c("white", "#80B1D3"))

geo <- Heatmap(mat,rect_gp = gpar(type = "none"),
               width = nrow(mat)*unit(6, "mm"),
               height = nrow(mat)*unit(6, "mm"),
               col=palette,na_col="white",
               bottom_annotation = c(bottom_group_ann),
               row_names_gp = gpar(fontsize = 8, fontface="italic"),
               column_names_gp = gpar(fontsize = 8),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name="Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                           labels_gp = gpar(fontsize = 8)),
               # cluster_rows = TRUE,
               # cluster_columns = TRUE,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(i >= j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                   grid.text(sprintf("%.f", mat[,1:nrow(mat)][i, j]), x, y, gp = gpar(fontsize = 6))
                 }
               }
)

# make fst heatmap
gene <- Heatmap(as.matrix(mat2[,1:nrow(mat2)]), rect_gp = gpar(type = "none"),
                width = nrow(mat2)*unit(6, "mm"),
                height = nrow(mat2)*unit(6, "mm"),
                right_annotation = row_group_ann,
                col=gene_col,na_col="grey",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name="FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                            labels_gp = gpar(fontsize = 8)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat2[,1:nrow(mat2)][i, j]), x, y, gp = gpar(fontsize = 6))
                  }
                }
                )

gene_width <- nrow(mat2)*unit(6, "mm")

draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Set the file name and parameters
filename <- "LantCama/outputs/LantCama_fst_plot.pdf"
width <- ncol(mat) * 0.28
height <- ncol(mat) * 0.28
dpi <- 300
units <- "cm"

# Set up the PNG device
pdf(filename, width = width, height = height)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Turn off the PNG device
dev.off()

####

ind_ho <- fastDiversity::individual_Ho(dms$gt, genetic_group_var=NULL, maf=0.02, max_missingness = 1)
m2$Ho <- ind_ho[match(m2$sample, names(ind_ho))]
ho_plot <- ggplot(m2, aes(x=cluster, y=Ho, fill=cluster))+
  geom_boxplot()+
  scale_fill_manual(values=tsne_cols)+
  theme_few()+
  theme(legend.position = 'none')+
  labs(x="HDBSCAN cluster", y="Ho (MAF 2%)")

fst_ho_combined <- ggarrange(fstp1+theme(legend.position = 'bottom'), ho_plot, widths=c(1,0.8),
                             # align='h',
                             labels="AUTO", font.label = list(face = "plain"))

fst_ho_combined

ggsave('LantCama/outputs/Figure2_fst_ho.png',dpi = 600, fst_ho_combined, width = 20, height = 8, units = "cm")


####

stat_site_cluster_maf0 <- fastDiversity::faststats(dms$gt %>% as.matrix(),
                                              genetic_group_variable = rep("Lcamara", length(dms$sample_names)),
                                              site_variable = dms$meta$site_cluster, 
                                              maf=0, max_missingness = 1)

stat_site_cluster_maf0 <- stat_site_cluster_maf0 %>%
  rename_with(~ ifelse(.x == "site", .x, paste0(.x, "_maf_0")))

stat_site_cluster_maf2 <- fastDiversity::faststats(dms$gt %>% as.matrix(),
                                                   genetic_group_variable = rep("Lcamara", length(dms$sample_names)),
                                                   site_variable = dms$meta$site_cluster, 
                                                   maf=0.02, max_missingness = 1)

stat_site_cluster_maf2 <- stat_site_cluster_maf2 %>%
  rename_with(~ ifelse(.x == "site", .x, paste0(.x, "_maf_2")))

library(dplyr)

summary_m2 <- m2 %>%
  group_by(site_cluster) %>%
  reframe(
    site = site,
    site_long = paste(unique(site_original), collapse = ":"), # Colon-delimited site_original values
    mean_lat = mean(lat, na.rm = TRUE),                      # Calculate mean latitude
    mean_long = mean(long, na.rm = TRUE),                    # Calculate mean longitude
    unique_clusters = paste(unique(cluster), collapse = ":"), # Colon-delimited cluster values
    unique_morphid2 = paste(unique(morphid2), collapse = ":") # Colon-delimited morphid2 values
  ) %>%
  unique()


stats_combined1 <- merge(summary_m2, stat_site_cluster_maf0, by.x="site_cluster", by.y="site")
stats_combined2 <- merge(stats_combined1, stat_site_cluster_maf2, stat_site_cluster_maf0, by.x="site_cluster", by.y="site")
cols_keep_stats <- c("site_cluster", "site","site_long", "mean_lat", "mean_long", "unique_clusters", 
                     "unique_morphid2", "n_maf_2", "Ho_maf_0", "uHe_maf_0", "uFis_maf_0", "loci_maf_0", 
                     "Ho_maf_2", "uHe_maf_2", "uFis_maf_2", "loci_maf_2"
)
stats_combined3 <- stats_combined2[,cols_keep_stats]
write.xlsx(stats_combined3,'LantCama/outputs/LantCama_supp_table_divstats.xlsx')
