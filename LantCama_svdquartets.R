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
# qc1       <- report_dart_qc_stats(d1, basedir, species, dataset, threshold_missing_loci = 0.8)

# d2        <- exclude.samples(d1,by="file", excluded_sample_file = "LantCama/meta/Lcam4_sampfilt_missing80.txt")

meta      <- read.meta.data.full.analyses.df(d1, basedir, species, dataset)
d3        <- dart.meta.data.merge(d1, meta) %>% remove.by.list(.,meta$sample_names[which(!is.na(meta$analyses[,'EA_AM']))])
d4        <- remove.poor.quality.snps(d3, min_repro=0.99, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
d6 <- remove.by.maf(d5, 0.02)
dms <- remove.by.missingness(d6, 0.8)

samples_to_keep_80 <- dms$sample_names
length(samples_to_keep_80)
loci_to_keep_80 <- colnames(dms$gt)
length(loci_to_keep_80)

write.table(samples_to_keep_80,'LantCama/meta/samples_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)
write.table(loci_to_keep_80,'LantCama/meta/loci_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)

# 
# meta      <- read.meta.data.full.analyses.df(dms, basedir, species, dataset)
# dms        <- dart.meta.data.merge(dms, meta) 

m2 <- d3$meta$analyses %>% as.data.frame

#### SVDq####

# dms_maf2 <- remove.by.maf(dms, 0.02)
# dms_maf5 <- remove.by.maf(dms, 0.05)
# 
# dms_1000 <- remove.loci.randomly(dms,1000)
# dms_5000 <- remove.loci.randomly(dms,5000)
# dms_10000 <- remove.loci.randomly(dms,10000)
# 
# dms_maf22 <- remove.by.missingness(dms_maf2, 0.8)

# dart2svdquartets(dms_maf2, RandRbase, species, dataset, add_pop=TRUE, pop=dms$sample_names)


svdq_pop_eilish <- dms$meta$analyses[,'svdq_pop']
svdq_pop_eilish[which(is.na(svdq_pop_eilish))] <- dms$sample_names[which(is.na(svdq_pop_eilish))]
dart2svdquartets(dms, RandRbase, species, dataset, add_pop=TRUE, pop=svdq_pop_eilish)
# dart2snapp(dms, RandRbase, species, dataset, add_pop=TRUE, pop=svdq_pop_eilish)

svdq_pop_pat <- dms$meta$analyses[,'svdq_pop2']
svdq_pop_pat[which(is.na(svdq_pop_pat))] <- dms$sample_names[which(is.na(svdq_pop_pat))]
dart2svdquartets(dms, RandRbase, species, dataset, add_pop=TRUE, pop=svdq_pop_pat)


#### iqtree ####

gl_unfiltered <- gl.read.dart(filename = "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/Report_DLan23-8067_SNP_mapping_2.csv", 
                              ind.metafile = '/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/LantCama_DLan23-8067_meta.csv')

samples_to_keep_80 <- read.csv('LantCama/meta/samples_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()
loci_to_keep_80 <- read.table('LantCama/meta/loci_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()

gl <- gl.keep.ind(gl_unfiltered,samples_to_keep_80) #keep only the listed individuals
indices <- which(gl$other$loc.metrics$AlleleID %in% loci_to_keep_80)
gl <- gl.keep.loc(gl,gl$loc.names[indices]) #keep only the listed lo

gl@pop <- gl@other[["ind.metrics"]][["svdq_pop_eilish"]]

glb <-   gl.read.dart(filename = "../BossFrag/BossFrag/dart_raw/Report-DBoss22-6915/Report_DBoss22-6915_1_moreOrders_SNP_2.csv")
gl2fasta(gl,outfile = "LantCama_eilish_iqtree.fasta", outpath = "LantCama/popgen",method =1,verbose = NULL)

##### SPLITSTREE ####

library(tanggle)
library(RSplitsTree)
# dms_1000 <- remove.by.missingness(dms_1000,0.8)
# splitstree(dist(dms$gt, method = "euclidean"), 'LantCama/outputs/all_LantCama_nexus_file_for_R_80missing_maf2.nex')
library(ggplot2)
library(phangorn)
library(ggforce)

# Read network data from Nexus file
Nnet <- phangorn::read.nexus.networx('LantCama/outputs/all_LantCama_nexus_file_for_R_80missing_maf2.nex')


x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)
x$svdq_pop[!is.na(x$sample)&is.na(x$svdq_pop)] <- "ungrouped"

net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)


levels2 <- unique(dms$meta$analyses[,"national2"]) %>% sort()
levels2_shapes <- setNames(1:nlevels(factor(levels2)), levels2)

hull <- x %>% group_by(svdq_pop) %>% 
  slice(chull(x, y))

hull <- hull[!hull$svdq_pop=='ungrouped',]

svdq_pop_colours <- named_list_maker(m2$svdq_pop, 'Spectral', 11)
svdq_pop_colours <- c(svdq_pop_colours, 'ungrouped'='grey30')
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")

x_filtered <- x[!is.na(x$national2),]

splitstree_plot_svdq <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.7, expand = 0.01, radius = 0.01,
             aes(fill = svdq_pop, color = "transparent")) +
  geom_point(data = x_filtered, aes(shape = national2), color="white", size=2) +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_filtered, aes(color = morphid2, shape = national2)) +
  scale_fill_manual(values = svdq_pop_colours, na.translate = FALSE) +
  scale_colour_manual(values = morphid_colours, na.translate = FALSE) +
  theme_void() +
  expand_limits(x = c(min(x_filtered$x) - 0.01 * net_x_axis, max(x_filtered$x) + 0.01 * net_x_axis),
                y = c(min(x_filtered$y) - 0.01 * net_y_axis, max(x_filtered$y) + 0.01 * net_y_axis)) +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  coord_fixed() +
  labs(color = "Morphotype", shape = "Origin", fill = "SVDq clusters") +
  guides(colour = guide_legend(title.position = "top", nrow = 2, override.aes = list(fill = NA, linetype = 0)), 
         fill = guide_legend(title.position = "top", nrow = 2, override.aes = list(color = NA)),
         shape = guide_legend(title.position = "top", nrow = 2))
# Save the plot
ggsave("LantCama/outputs/LantCama_splitstree_svdqpops_80miss_maf2.pdf",
       splitstree_plot_svdq, width = 20, height = 20, units = "cm", dpi = 600)

ggsave("LantCama/outputs/LantCama_splitstree_svdqpops_80miss_maf2.png",
       splitstree_plot_svdq, width = 20, height = 20, units = "cm", dpi = 600)

###### TREE ########

library(phangorn)

#hamming distance
genetic_distances <- dist(dms$gt, method = "binary")

bootstrap_trees <- bootstrap.phyDat(genetic_distances, FUN = function(x) upgma(dist.dna(x)))
# genetic_distances <- dist(dms$gt, method = "euclidean")

# genetic_distances <- ape::nei.dist(as.matrix(dms$gt))

tree <- upgma(genetic_distances)
# tree <- nj(genetic_distances)
# tree <- wpgma(genetic_distances)

ggtree_obj <- ggtree(tree)

x1 <- m2[m2$sample %in% ggtree_obj$data$label,]
rownames(x1) <- x1$sample


ggtree_obj <- ggtree(tree, size=0.01) %<+% x1 
# geom_tippoint(aes(color = cluster, shape=national2), size=0.7) + theme_tree2()
ggtree_obj

cc <-   named_list_maker(x1$cluster_histogram, 'Spectral',11)
cc2 <- named_list_maker(x1$cluster, 'Paired',11)
cc2 <- named_list_maker(x1$national2, 'Paired',11)

hmt <-  gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop2', 'svdq_pop','national2')]),
                 offset=0.005, width=.1,font.size=2,
                 colnames_angle=90, colnames_position="top")+
  # scale_fill_manual(values=named_list_maker(x$svdq_pop, 'Paired',11))+
  theme_tree2()+
  # theme(legend.position = "right")+
  scale_y_continuous(expand=c(0.02, 0))+
  geom_vline(xintercept = 0.05, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0.1, color="grey80", linetype="dotted")+
  geom_rootedge(0.001, size=0.01)+
  geom_tiplab(aes(label = label), size=0.75)  # Add this line to include tip labels
ggsave("LantCama/outputs/LantCama_euclidean_upgma_maf2.pdf",
       hmt, width = 20, height = 30, units = "cm", dpi=600)


hmt <-  gheatmap(ggtree_obj, x1[,c('cluster_histogram','morphid2','national2')],
                 offset=0.005, width=.1,font.size=2,colnames_offset_y=c(10,15,10),
                 custom_column_labels=c("Cluster","Morphotype","Country"),
                 colnames_angle=90, colnames_position="top")+
  scale_fill_manual(values=c(cc,cc2, morphid_colours))+
  theme_tree2()+
  theme(legend.position = "right")+
  scale_y_continuous(expand=c(0.02, 0))+
  geom_vline(xintercept = 0.05, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0.1, color="grey80", linetype="dotted")+
  geom_rootedge(0.001, size=0.01)+
  geom_tiplab(aes(label = label), size=0.75)  # Add this line to include tip labels

ggsave("LantCama/outputs/LantCama_euclidean_upgma_all_maf2.pdf",
       hmt, width = 30, height = 35, units = "cm", dpi=600)

#### pca ####

gen_d5 <- new("genlight", dms[["gt"]]) #convert df to genlight object for glPca function
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
  guides(colour = guide_legend(title.position = "top"))+#+
  scale_colour_manual(values=morphid_colours)
# scale_shape_manual(values=cluster2ecies_shapes)
pca_plot1

pca_plot2 <- ggplot(g_pca_df2, aes(x=PC3, y=PC4, colour=cluster2, shape=national2))+ xlab(pcnames[3])+ylab(pcnames[4])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))
# scale_colour_manual(values=cluster2ecies_colours)+
# scale_shape_manual(values=cluster2ecies_shapes)

pca_plot3 <- ggplot(g_pca_df2, aes(x=PC5, y=PC6, colour=cluster2, shape=national2))+ xlab(pcnames[5])+ylab(pcnames[6])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))
# scale_colour_manual(values=cluster2ecies_colours)+
# scale_shape_manual(values=cluster2ecies_shapes)

all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3, labels=c("A","B","C"),
                            common.legend = TRUE, ncol=3, legend = "bottom")
all3_pca_plots



# ###
# library(dartR)
# gl_unfiltered <- gl.read.dart(filename = "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/Report_DLan23-8067_SNP_mapping_2.csv")
# 
# samples_to_keep_80 <- read.csv('LantCama/meta/samples_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()
# loci_to_keep_80 <- read.table('LantCama/meta/loci_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()
# 
# # gl <- gl.drop.pop() #remove listed populations
# # gl <- gl.keep.pop() #keep only the listed populations
# # gl <- gl.drop.ind() #remove listed individuals
# gl <- gl.keep.ind(gl_unfiltered,samples_to_keep_80) #keep only the listed individuals
# # gl <- gl.drop.loc() #remove listed loci
# indices <- which(gl$other$loc.metrics$AlleleID %in% loci_to_keep_80)
# 
# gl <- gl.keep.loc(gl,gl$loc.names[indices]) #keep only the listed lo
# gl$other$loc.metrics$TrimmedSequence <- gl$other$loc.metrics$TrimmedSequenceSnp[indices]
# gl2fasta(gl,outfile = "output.fasta",outpath='/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen', method=3)
# 
# sild <- read.csv('/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/Report-DLan23-8067/Report_DLan23-8067_4_moreOrders_SilicoDArT_1.csv', skip=6)
# View(sild)
# sild <- sild[sild$Reproducibility>0.96,]
# nrow(sild)
# 
# sild2 <- sild[,11:ncol(sild)]
# sild2[sild2=='-'] <- NA
# sild2 <- as.matrix(sild2)
# unique(sild2)
# rowMeans(sild2, na.rm=TRUE)
