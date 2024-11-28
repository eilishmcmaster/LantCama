library(ggsignif) 
library(ggforce)
library(RRtools)
library(ade4)
library(adegenet)
library(ape)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(diveRsity)
library(dplyr)
library(geosphere)
library(ggfortify)
library(ggmap)
library(ggrepel)
library(fastDiversity)
library(ggpubr)
library(ggthemes)
library(ggtree)
library(heatmaply)
library(lattice)
library(openxlsx)
library(ozmaps)
library(RColorBrewer)
library(RRtools)
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
d4        <- remove.poor.quality.snps(d3, min_repro=0.99, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
dms <- remove.by.missingness(d5, 0.8)

m2 <- d3$meta$analyses %>% as.data.frame

#### colours 
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")
svdq_pop_colours <- named_list_maker(m2$svdq_pop, 'Spectral', 11)
svdq_pop_colours <- c(svdq_pop_colours, 'ungrouped'='grey30')


#### pca ####
dms_maf2 <- remove.by.maf(dms, 0.02)

gen_d5 <- new("genlight", dms_maf2[["gt"]]) #convert df to genlight object for glPca function
gen_pca <- glPca(gen_d5, parallel=TRUE, nf=6) #do pca -- this method somehow allows the input to hav1 NAs

g_pca_df <- gen_pca[["scores"]] #extract PCs
g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # some in DArT are not in meta?

pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:6]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes

pca_plot1 <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=morphid2))+ xlab(pcnames[1])+ylab(pcnames[2])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
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

####
library(dbscan)

eps <- 14
min_pts <- 10

genotype_matrix <- dms$gt %>% as.matrix()
d <- dist(genotype_matrix, method = "euclidean")

d_matrix <- as.matrix(d)  # Convert dist object to matrix

# Extract the upper triangle
upper_triangle <- d_matrix[upper.tri(d_matrix)]

# Plot the histogram
hist(upper_triangle, 
     breaks = 100,  # Adjust number of bins as needed
     main = "Histogram of Upper Triangle of Distance Matrix", 
     xlab = "Pairwise Distance", 
     col = "skyblue", 
     border = "white")
abline(v = eps, col = "red", lty = 2, lwd = 2)  # lty = 2 for dotted line, lwd = 2 for thickness

kNNdistplot(d, k = min_pts) # k = min pts
abline(h = eps, col = "red", lty = 2)  # Adjust h to where the elbow appears

# Use DBSCAN to cluster with a strict distance constraint
# db_result <- dbscan(d, eps = eps, minPts = min_pts)  
db_result <- hdbscan(d, minPts = min_pts)  

print(db_result$cluster)
clusters <- db_result$cluster

db_df <- data.frame(sample=names(d),db_cluster=clusters)
db_df$db_cluster[which(db_df$db_cluster==0)] <- NA

##

library(ComplexHeatmap)
library(circlize)

# Create a color mapping for clusters
cluster_colors <- structure(
  c(rainbow(max(db_result$cluster, na.rm = TRUE))), # Gray for noise
  names = c(seq_len(max(db_result$cluster, na.rm = TRUE)))
)

# Create row and column annotations
row_anno <- rowAnnotation(
  db_cluster = as.factor(db_df$db_cluster),
  col = list(db_cluster = cluster_colors),
  annotation_legend_param = list(
    db_cluster = list(title = "DBSCAN Cluster")
  )
)

col_anno <- HeatmapAnnotation(
  db_cluster = as.factor(db_df$db_cluster),
  col = list(db_cluster = cluster_colors),
  annotation_legend_param = list(
    db_cluster = list(title = "DBSCAN Cluster")
  ),
  which = "column"
)

# Create the heatmap
ht <- Heatmap(
  d_matrix,
  name = "Distance",
  clustering_method_columns = "average",
  clustering_method_rows = "average",
  top_annotation = col_anno,
  right_annotation = row_anno,
  
  col = colorRamp2(c(min(d_matrix), max(d_matrix)), c("white", "blue"))
)
draw(ht, merge_legends = TRUE)

### TSNE ####

library(Rtsne)

library(dbscan)
library(Rtsne)  # For t-SNE
library(ggplot2)  # For plotting

# Existing code to compute the distance matrix
genotype_matrix <- dms$gt %>% as.matrix()
d <- dist(genotype_matrix, method = "euclidean")
d_matrix <- as.matrix(d)

# Perform t-SNE on the distance matrix (or you could use the original genotype matrix)
tsne_result <- Rtsne(d_matrix, dims = 2, perplexity = 10, 
                     check_duplicates = FALSE, theta=0,
                     is_distance=TRUE,
                     num_threads=4)

# Run HDBSCAN on the t-SNE output (2D coordinates)
tsne_df <- data.frame(tsne_result$Y)  # Extract the 2D coordinates from t-SNE
colnames(tsne_df) <- c("tSNE1", "tSNE2")  # Rename columns for clarity

tsne_dist <- dist(tsne_result$Y) %>% as.matrix()

# Extract the upper triangle
upper_triangle <- tsne_dist[upper.tri(tsne_dist)]

# Plot the histogram
hist(upper_triangle, 
     breaks = 100,  # Adjust number of bins as needed
     main = "Histogram of Upper Triangle of Distance Matrix", 
     xlab = "Pairwise Distance", 
     col = "skyblue", 
     border = "white")
abline(v = 10, col = "red", lty = 2, lwd = 2)  # lty = 2 for dotted line, lwd = 2 for thickness


# Apply HDBSCAN to the t-SNE coordinates (2D space)
hdbscan_result <- hdbscan(tsne_df, minPts = 5)
# hdbscan_result <- dbscan(tsne_df, minPts = 10, eps=4)

hdb_scan_qc <- data.frame(cluster=hdbscan_result$cluster, prob=hdbscan_result$membership_prob,
                   cd = hdbscan_result$coredist, os = hdbscan_result$outlier_scores)

hdb_scan_qc_summary_stats <- hdb_scan_qc %>%
  # filter(cluster != 0) %>%  # Exclude cluster 0
  group_by(cluster) %>%
  summarise(
    mean_os = mean(os, na.rm = TRUE),  # Compute mean of outlier scores, ignoring NA values
    sd_os = sd(os, na.rm = TRUE)       # Compute standard deviation of outlier scores
  ) %>%
  mutate( # Z-score indicates how many standard deviations an element is from the mean. 
    z_score = (mean_os - mean(mean_os)) / sd(mean_os)
  )

ggplot()+
  geom_histogram(data=hdb_scan_qc, mapping=aes(x=cd))+
  facet_wrap(~cluster)

ggplot()+
  geom_histogram(data=hdb_scan_qc, mapping=aes(x=os))+
  facet_wrap(~cluster)+
  geom_vline(data = hdb_scan_qc_summary_stats, aes(xintercept = mean_os+1.5*(sd_os)), linetype = "dashed", color = "red")

hdb_scan_qc_summary_stats

unusually_high_clusters <- hdb_scan_qc_summary_stats %>%
  filter(z_score > 1.5)

hdbscan_result$cluster[which(hdbscan_result$cluster %in% unusually_high_clusters$cluster)] <- 0 # remove unusual clusters

# Add the cluster labels to the t-SNE data
tsne_df$db_cluster <- hdbscan_result$cluster
tsne_df$db_cluster[tsne_df$db_cluster == 0] <- NA  # Assign NA to noise points

tsne_cols <- c(rainbow(length(unique(tsne_df$db_cluster))))
names(tsne_cols) <- unique(tsne_df$db_cluster)
# Plot the t-SNE result with HDBSCAN clusters

hull <- tsne_df %>% 
  filter(!is.na(db_cluster)) %>%
  group_by(db_cluster) %>% 
  slice(chull(tSNE1, tSNE2))

ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = factor(db_cluster))) +
  geom_shape(data=hull, mapping=aes(group=db_cluster), color="transparent", alpha=0.1,
             expand = 0.01, radius=0.01)+
  geom_point(alpha=0.5) +
  scale_color_manual(values = tsne_cols) +
  labs(title = "t-SNE with HDBSCAN Clusters",
       x = "t-SNE 1", y = "t-SNE 2", color = "Cluster") +
  theme_bw()


#
library(ComplexHeatmap)
library(circlize)

svdq_pop_colours <- colorRampPalette(c("grey20", "grey90"))(7)
names(svdq_pop_colours) <- 1:7
svdq_pop_colours <- c(c(eacp="#AA3377", per1="#228833",  eawt="#66CCEE"),svdq_pop_colours)

clusters2 <- hdbscan_result$cluster

db_df2 <- data.frame(sample=names(d),db_cluster=clusters2)
db_df2$db_cluster[which(db_df2$db_cluster==0)] <- NA
db_df2 <- merge(db_df2, m2, by='sample')
##

d_matrix2 <- merge(d_matrix, db_df2, by.x=0, by.y='sample')
rownames(d_matrix2) <- d_matrix2$Row.names
d_matrix2$Row.names <- NULL
d_matrix2 <- d_matrix2[match(colnames(d_matrix2)[1:nrow(d_matrix2)],rownames(d_matrix2)),]

tsne_cols2 <- tsne_cols[!is.na(names(tsne_cols))]

# Create row and column annotations
row_anno <- rowAnnotation(
  db_cluster = as.factor(d_matrix2$db_cluster),
  col = list(db_cluster = tsne_cols2),
  annotation_legend_param = list(
    db_cluster = list(title = "DBSCAN Cluster")
  ),
  annotation_name_gp = gpar(fontsize = 0)
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
  svdq_pop = d_matrix2$svdq_pop_label,
  col = list(svdq_pop = svdq_pop_colours),
  annotation_legend_param = list(
    svdq_pop = list(title = "SVDQ pop label")
  ),
  annotation_name_gp = gpar(fontsize = 0)
)

col_anno <- HeatmapAnnotation(
  db_cluster = as.factor(d_matrix2$db_cluster),
  col = list(db_cluster = tsne_cols2),
  annotation_legend_param = list(
    db_cluster = list(title = "DBSCAN Cluster")
  ),
  which = "column"
)

# Create the heatmap
ht <- Heatmap(
  d_matrix2[,1:nrow(d_matrix2)] %>% as.matrix,
  name = "Distance",
  clustering_method_columns = "average",
  clustering_method_rows = "average",
  top_annotation = col_anno,
  right_annotation = c(row_anno, row_anno2,row_anno3),
  show_column_names = FALSE,  # Remove column names
  show_row_names = FALSE,   
  col = colorRamp2(c(min(d_matrix), max(d_matrix)), c("white", "blue"))
)
draw(ht, merge_legends = TRUE)

