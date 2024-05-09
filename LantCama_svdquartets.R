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

d2        <- exclude.samples(d1,by="file",
                             excluded_sample_file = "LantCama/meta/Lcam4_sampfilt_missing80.txt", remove_fixed_loci = TRUE)
d2        <- remove.poor.quality.snps(d2, min_repro=0.96, max_missing=0.8)
# qc2       <- report_dart_qc_stats(d2, basedir, species, dataset)
# d2 <- remove.by.missingness(d2, 0.5)
d3        <- sample.one.snp.per.locus.random(d2, seed=12345) 
# qc3       <- report_dart_qc_stats(d3, basedir, species, dataset)

v         <- "RRv0003"

meta      <- read.meta.data.full.analyses.df(d3, basedir, species, dataset)
# meta      <- read_meta_info(d3, basedir, species, dataset, version=v) 

dms        <- dart.meta.data.merge(d3, meta)
# dm        <- merge_gt_meta_data(d3, meta)
m2 <- dms$meta$analyses %>% as.data.frame
# dmv       <- arrange_data_by_analysis_field(dm, "EA_only", basedir, species, dataset)


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

ggsave("LantCama/outputs/LantCama_splitstree_nation_morpho.svg",
       splitstree_plot, width = 20, height = 30, units = "cm", dpi=600)

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

ggsave("LantCama/outputs/LantCama_splitstree_cluster.svg",
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

