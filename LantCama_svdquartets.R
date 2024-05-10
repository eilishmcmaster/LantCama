library(RRtools)
source('https://github.com/eilishmcmaster/SoS_functions/blob/33bd7065baa91b7e2a1800b6481d63573fb38d28/dart2svdquartets.r?raw=TRUE')

topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory 
missingness <- 0.3
species   <- "LantCama"
dataset   <- "DLan23-8067"

basedir <- ""

d1        <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)
# qc1       <- report_dart_qc_stats(d1, basedir, species, dataset, threshold_missing_loci = 0.8)

d2        <- exclude.samples(d1,by="file", excluded_sample_file = "LantCama/meta/Lcam4_sampfilt_missing80.txt")

meta      <- read.meta.data.full.analyses.df(d2, basedir, species, dataset)
d3        <- dart.meta.data.merge(d2, meta) %>% remove.by.list(.,meta$sample_names[which(!is.na(meta$analyses[,'EA_AM']))])
d3        <- remove.poor.quality.snps(d3, min_repro=0.96, max_missing=0.4)%>% remove.fixed.snps()
dms        <- sample.one.snp.per.locus.random(d3, seed=12345) 
dmsx <- remove.by.missingness(dms, 0.5)
samples_to_keep_80 <- dms$sample_names
loci_to_keep_80 <- colnames(dms$gt)

write.table(samples_to_keep_80,'LantCama/meta/samples_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)
write.table(loci_to_keep_80,'LantCama/meta/loci_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)

m2 <- dms$meta$analyses %>% as.data.frame

#### SVDq####
# dart2svdquartets(dms, RandRbase, species, dataset, add_pop=TRUE, pop=dms$sample_names)


####

#### hclust 
# 
# hclust_avg <- hclust(dist(dms$gt), method = 'average')
# plot(hclust_avg)
# library(ggtree)
library(pvclust)#https://cran.r-project.org/web/packages/pvclust/index.html
# 
# result <- pvclust(dms$gt, method.dist="euclidean", method.hclust="average", nboot=10, parallel=TRUE)


##### SPLITSTREE ####

library(tanggle)
library(RSplitsTree)

# splitstree(dist(dms$gt), 'LantCama/outputs/all_LantCama_nexus_file_for_R_80%missing.nex')

#need to open and save the file in Splitstree app for it to open here, IDK why
Nnet <- phangorn::read.nexus.networx('LantCama/outputs/all_LantCama_nexus_file_for_R_80%missing.nex')

x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)

net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)



# levels <- unique(dms$meta$analyses[,"morphid"]) %>% sort()
# # Get the default ggplot colors for the levels
# colors <- scales::hue_pal()(length(levels))
# # Create a named list of colors
# species_colours <- setNames(colors, levels)

levels2 <- unique(dms$meta$analyses[,"national2"]) %>% sort()
levels2_shapes <- setNames(1:nlevels(factor(levels2)), levels2)


#set up colours (Tol colour palette)
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")


splitstree_plot <- ggplot(Nnet, mapping = aes_(~x, ~y), layout = "slanted", mrsd = NULL, 
                          as.Date = FALSE, yscale = "none", yscale_mapping = NULL, 
                          ladderize = FALSE, right = FALSE, branch.length = "branch.length", 
                          ndigits = NULL)+
  geom_splitnet(layout = "slanted", size=0.2, color='grey20')+
  geom_point(data=x, aes(x, y, colour=morphid2, shape=national2), size=1, fill="white")+#shape=national2
  scale_colour_manual(values=morphid_colours, na.translate=FALSE, name='Morphotype')+
  scale_shape_manual(values=c(16,17), na.translate=FALSE, name='Country')+
  # geom_tiplab2(aes(label=label), size=2, hjust=-0.2)+
  theme_void()+
  expand_limits(x=c(min(x$x)-0.01*net_x_axis, max(x$x)+0.01*net_x_axis),
                y=c(min(x$y)-0.01*net_y_axis, max(x$y)+0.01*net_y_axis))+
  theme(legend.position = "bottom", legend.key.size = unit(0.5, 'lines'))+coord_fixed()+
  guides(colour = guide_legend(title.position = "top", nrow = 3), shape = guide_legend(title.position = "top", nrow = 4))

ggsave("LantCama/outputs/LantCama_splitstree_nation_morpho.pdf",
       splitstree_plot, width = 20, height = 20, units = "cm", dpi=600)

splitstree_plot_cluster <- ggplot(Nnet, mapping = aes_(~x, ~y), layout = "slanted", mrsd = NULL, 
       as.Date = FALSE, yscale = "none", yscale_mapping = NULL, 
       ladderize = FALSE, right = FALSE, branch.length = "branch.length", 
       ndigits = NULL)+
  geom_splitnet(layout = "slanted", size=0.2)+
  geom_point(data=x, aes(x, y, colour=cluster_histogram))+#shape=national2
  scale_colour_manual(values= scales::hue_pal()(length(unique(x$cluster_histogram))), na.translate=FALSE, name='Morphotype')+
  # scale_shape_manual(values=levels2_shapes, na.translate=FALSE, name='Country')+
  # geom_tiplab2(aes(label=label), size=2, hjust=-0.2)+
  theme_void()+labs(color="", shape="")+
  expand_limits(x=c(min(x$x)-0.01*net_x_axis, max(x$x)+0.01*net_x_axis),
                y=c(min(x$y)-0.01*net_y_axis, max(x$y)+0.01*net_y_axis))+
  theme(legend.position = "bottom")+coord_fixed()+
  guides(colour = guide_legend(title.position = "top", nrow = 3), shape = guide_legend(title.position = "top", nrow = 4))

ggsave("LantCama/outputs/LantCama_splitstree_cluster.pdf",
       splitstree_plot_cluster, width = 20, height = 30, units = "cm", dpi=600)
# 
# 
# tree_tree_plot <- ggarrange(splitstree_plot, ggtree_plot, ncol=2, labels=c("D", "E"), legend="none")
# pca_tree_tree_plot <- ggarrange(all3_pca_plots,NULL, tree_tree_plot, nrow=3, heights=c(1.25,0.1,2))&
#   theme(plot.background = element_rect(fill = "white"))+border("white")
# pca_tree_tree_plot
# 
# 
# ggsave("LantCama/outputs/paper/fig2_LantCama_pca_tree_tree.pdf",
#        pca_tree_tree_plot, width = 18, height = 18, units = "cm", dpi=600)

###### TREE ########

library(phangorn)

#hamming distance
genetic_distances <- dist(dmsx$gt, method = "binary")
tree <- upgma(genetic_distances)
# tree <- nj(genetic_distances)
# tree <- wpgma(genetic_distances)

ggtree_obj <- ggtree(tree)

x1 <- m2[m2$sample %in% ggtree_obj$data$label,]
rownames(x1) <- x1$sample

# Assuming you already have the ggtree object ggtree_obj and tree object
# 
# # Plot the tree with modified branch width and distances on x-axis
# ggtree_obj <- ggtree(tree, layout="circular") %<+% x1 +
#   geom_tippoint(aes(color = morphid2, shape=national2), size=1) +  
#   # scale_colour_manual(values=morphid_colours, na.translate=FALSE, name='Morphotype')+
#   scale_shape_manual(values=c(16,17), na.translate=FALSE, name='Country')
# 
# ggsave("LantCama/outputs/LantCama_hamming_upgma.pdf",
#        ggtree_obj, width = 30, height = 60, units = "cm", dpi=600)

ggtree_obj <- ggtree(tree, size=0.01) %<+% x1 
# geom_tippoint(aes(color = cluster, shape=national2), size=0.7) + theme_tree2()
ggtree_obj

cc <-   named_list_maker(x1$cluster_histogram, 'Spectral',11)
cc2 <- named_list_maker(x1$cluster, 'Paired',11)
cc2 <- named_list_maker(x1$national2, 'Paired',11)

hmt <-  gheatmap(ggtree_obj, x1[,c('cluster_histogram','morphid2','national2')],
                 offset=0, width=.1,font.size=2,colnames_offset_y=c(10,15,10),
                 custom_column_labels=c("Cluster","Morphotype","Country"),
                 colnames_angle=90, colnames_position="top")+
  scale_fill_manual(values=c(cc,cc2, morphid_colours))+
  theme_tree2()+theme(legend.position = "none")+
  scale_y_continuous(expand=c(0.02, 0))+
  geom_vline(xintercept = 0.05, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0, color="grey80", linetype="dotted")+
  geom_vline(xintercept = 0.1, color="grey80", linetype="dotted")+
  geom_rootedge(0.001, size=0.01)

ggsave("LantCama/outputs/LantCama_hamming_upgma.pdf",
       hmt, width = 20, height = 30, units = "cm", dpi=600)


####

gl_unfiltered <- gl.read.dart(filename = "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/Report_DLan23-8067_SNP_mapping_2.csv")

samples_to_keep_80 <- read.csv('LantCama/meta/samples_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()
loci_to_keep_80 <- read.table('LantCama/meta/loci_to_keep_80%.csv', header=FALSE) %>% t() %>% as.vector()

# gl <- gl.drop.pop() #remove listed populations
# gl <- gl.keep.pop() #keep only the listed populations
# gl <- gl.drop.ind() #remove listed individuals
gl <- gl.keep.ind(gl_unfiltered,samples_to_keep_80) #keep only the listed individuals
# gl <- gl.drop.loc() #remove listed loci
indices <- which(gl$other$loc.metrics$AlleleID %in% loci_to_keep_80)

gl <- gl.keep.loc(gl,gl$loc.names[indices]) #keep only the listed lo
gl$other$loc.metrics$TrimmedSequence <- gl$other$loc.metrics$TrimmedSequenceSnp[indices]
gl2fasta(gl,outfile = "output.fasta",outpath='/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen', method=3)

sild <- read.csv('/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/Report-DLan23-8067/Report_DLan23-8067_4_moreOrders_SilicoDArT_1.csv', skip=6)
View(sild)
sild <- sild[sild$Reproducibility>0.96,]
nrow(sild)

sild2 <- sild[,11:ncol(sild)]
sild2[sild2=='-'] <- NA
sild2 <- as.matrix(sild2)
unique(sild2)
rowMeans(sild2, na.rm=TRUE)
