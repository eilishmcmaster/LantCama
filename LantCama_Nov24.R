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
genotype_matrix <- dms$gt %>% as.matrix()
d <- dist(genotype_matrix, method = "euclidean")
# hc <- hclust(d, method = "ward.D2")
# hc_upgma <- hclust(d, method = "average")
# plot(hc)  # Visualize dendrogram
# plot(hc_upgma)

# Ensure `d` is a distance object
d_matrix <- as.matrix(d)  # Convert dist object to matrix

# Extract the upper triangle
upper_triangle <- d_matrix[upper.tri(d_matrix)]

eps <- 10

# Plot the histogram
hist(upper_triangle, 
     breaks = 100,  # Adjust number of bins as needed
     main = "Histogram of Upper Triangle of Distance Matrix", 
     xlab = "Pairwise Distance", 
     col = "skyblue", 
     border = "white")
abline(v = eps, col = "red", lty = 2, lwd = 2)  # lty = 2 for dotted line, lwd = 2 for thickness


kNNdistplot(d, k = 5) # k = min pts
abline(h = eps, col = "red", lty = 2)  # Adjust h to where the elbow appears
db_result <- dbscan(d, eps = eps, minPts = 5)
print(db_result)

db_df <- data.frame(sample=names(d),db_cluster=db_result$cluster)
db_df$db_cluster[which(db_df$db_cluster==0)] <- NA

g_pca_df3 <- merge(g_pca_df2,db_df, by.x='Row.names',by.y='sample')

pca_plot1 <- ggplot(g_pca_df3, aes(x=PC1, y=PC2, colour=factor(db_cluster)))+ xlab(pcnames[1])+ylab(pcnames[2])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))
  # scale_colour_manual(values=morphid_colours)
pca_plot1



# 
# set.seed(123)  # For reproducibility
# kmeans_result <- kmeans(d, centers = 3)
# clusters <- kmeans_result$cluster
# pam_result <- pam(d, k = 3)  # Choose the number of clusters
# clusters <- pam_result$clustering
# heatmap(d, Rowv = as.dendrogram(hc), symm = TRUE)
# 
# pca <- prcomp(genotype_matrix, scale. = TRUE)
# plot(pca$x[, 1:2], col = clusters, pch = 16, main = "PCA Clustering")
# 
# sil <- silhouette(clusters, dist = d)
# plot(sil)

# UMAP ###################################################################
# UMAP reduces dimensions by retaining relationships between points (sample focussed)
# UMAP calculations
library(umap)
custom.config = umap.defaults
custom.config$random_state = 330

umer <- umap(d %>% as.matrix(), config=custom.config, n_neighbors = 5, min_dist = 0.05) # run umap

umap_df <- umer$layout %>% as.data.frame() #extract output vectors
umap_df2 <- merge(umap_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=TRUE) #add metadata
umap_df2 <- merge(umap_df2, db_df, by.x='Row.names',by.y='sample')

# hull2 <- umap_df2 %>% group_by(svdq_pop_label) %>% 
#   slice(chull(V1, V2))
# hull2 <- hull2[!hull2$svdq_pop_label=='ungrouped',]
hull2 <- umap_df2 %>% group_by(db_cluster) %>% 
  slice(chull(V1, V2))
hull2 <- hull2[hull2$db_cluster!=0,]
# 
# # Find the point with the largest V1 for each svdq_pop_label
# label_points2 <- umap_df2 %>%
#   group_by(svdq_pop_label) %>%
#   filter(V1 == max(V1, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Filter out the "ungrouped" label if needed
# label_points2 <- label_points2 %>% filter(svdq_pop_label != 'ungrouped')

# Create the plot
umap_plot <- ggplot(umap_df2, aes(x = V1, y = V2, colour = morphid2)) +
  geom_shape(data = hull2, alpha = 0.5, expand = 0.01, radius = 0.01,
             aes(group = db_cluster), color = "transparent", show.legend=FALSE) +
  # scale_fill_manual(values = svdq_pop_colours, na.translate = FALSE, na.value = "transparent") +
  geom_point() +
  # geom_text_repel(
  #   data = label_points2, aes(x = V1, y = V2, label = svdq_pop_label),
  #   size = 3,  color = "black", nudge_x = 0.3, nudge_y=0.4
  # ) +
  theme_bw() +
  scale_colour_manual(values = morphid_colours) +
  guides(colour = guide_legend(
    title.position = "top",
    override.aes = list(fill = NA, linetype = 0)
  )) +
  theme(
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(colour = "Morphotype", fill = "Clusters", x = "", y = "")

# Display the plot
print(umap_plot)

