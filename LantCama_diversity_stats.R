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
# qc1       <- report_dart_qc_stats(d1, basedir, species, dataset, threshold_missing_loci = 0.8)

# d2        <- exclude.samples(d1,by="file", excluded_sample_file = "LantCama/meta/Lcam4_sampfilt_missing80.txt")

meta      <- read.meta.data.full.analyses.df(d1, basedir, species, dataset)
d3        <- dart.meta.data.merge(d1, meta) %>% remove.by.list(.,meta$sample_names[which(!is.na(meta$analyses[,'EA_AM']))])
d4        <- remove.poor.quality.snps(d3, min_repro=0.99, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12342) 

samples_to_keep_80 <- read.csv('LantCama/meta/samples_to_keep_80%.csv', header = FALSE) %>% unlist() %>% as.vector()

dms <- remove.by.list(d5, samples_to_keep_80)

d6 <- remove.by.maf(d5, 0.02)
dms_maf2 <- remove.by.missingness(d6, 0.8)

m2 <- d3$meta$analyses %>% as.data.frame

svdq_pop_colours <- named_list_maker(m2$svdq_pop, 'Spectral', 11)
svdq_pop_colours <- c(svdq_pop_colours, 'ungrouped'='grey30')
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")
nation_colours <- named_list_maker(m2$national2, 'Paired',11)


##### modified function, more dps 

faststats <- function (gt, genetic_group_variable, site_variable, minimum_n = 3, 
                       minimum_loci = 50, maf = 0.05, max_missingness = 0.3) 
{  out_matrix <- matrix(NA, 0, 10)
  colnames(out_matrix) <- c("genetic_group", "site", "Ar", 
                            "Ho", "He", "uHe", "Fis", "uFis", "loci", "n")
  
  genetic_group_freq <- table(genetic_group_variable)
  genetic_groups <- names(which(genetic_group_freq >= minimum_n))
  
  for (group in genetic_groups) {
    gt_group <- gt[which(genetic_group_variable == group), 
    ]
    site_variable_group <- site_variable[which(genetic_group_variable == 
                                                 group)]
    not_missing_loci <- which(colMeans(is.na(gt_group)) <= 
                                max_missingness)
    gt_group_missing <- gt_group[, not_missing_loci]
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
    passing_maf_loci <- which(loci_mafs >= maf)
    gt_group_missing_maf <- gt_group_missing[, passing_maf_loci]
    
    if (length(passing_maf_loci) < minimum_loci) {
      print(paste(group, "does not have enough loci (",
                  length(passing_maf_loci), "loci: minimum is,",
                  minimum_loci, ")"))
      print("Proceeding to next genetic group")
      (next)()
    }
    site_freq <- table(site_variable_group)
    sites <- names(which(site_freq >= minimum_n))
    site_variable_group <- site_variable_group[which(site_variable_group %in% 
                                                       sites)]
    group_out_matrix <- matrix(NA, length(sites), 10)
    for (s in 1:length(sites)) {
      site <- sites[s]
      gt_site <- gt_group_missing_maf[which(site_variable_group == 
                                              site), ]
      Ho <- calculate_Ho(gt_site)
      Hes <- calculate_Hes(gt_site)
      He <- mean(Hes, na.rm = TRUE)
      ns <- colSums(!is.na(gt_site))
      uHe <- calculate_uHe(ns, Hes)
      Fis <- 1 - (Ho/He)
      uFis <- 1 - (Ho/uHe)
      loci <- ncol(gt_site)
      Ar <- calculate_Ar(gt_site)
      n <- nrow(gt_site)
      group_out_matrix[s, ] <- c(group, site, round(Ar, 
                                                    5), round(Ho, 5), round(He, 5), round(uHe, 5), 
                                 round(Fis, 5), round(uFis, 5), loci, n)
    }
    out_matrix <- rbind(out_matrix, group_out_matrix)
    print(paste(group, "complete"))
  }
  print("Process complete!")
  return(out_matrix %>% as.data.frame())
}


#### stats ####

pops_with_na <- dms$meta$analyses[,'svdq_pop_label']
pops_with_na[is.na(pops_with_na)] <- "Ungrouped"

s3 <- faststats(dms$gt, genetic_group_variable = rep('lantana camara', length(dms$sample_names)), site_variable = pops_with_na, minimum_n = 2, maf = 0.00, max_missingness = 1)

s4 <- faststats(dms_maf2$gt, genetic_group_variable = rep('lantana camara', length(dms_maf2$sample_names)), site_variable = pops_with_na, minimum_n = 2, maf = 0, max_missingness = 1)


# pop_morphEA columns
s5 <- faststats(dms$gt, genetic_group_variable = rep('lantana camara', length(dms$sample_names)), site_variable = dms$meta$analyses[,'pop_morphEA'], minimum_n = 2, maf = 0, max_missingness = 1)

s6 <- faststats(dms_maf2$gt, genetic_group_variable = rep('lantana camara', length(dms_maf2$sample_names)), site_variable = dms_maf2$meta$analyses[,'pop_morphEA'], minimum_n = 2, maf = 0, max_missingness = 1)


summary_df <- dms$meta$analyses %>% as.data.frame() %>%
  group_by(morphid2, pop_morphEA) %>%#svdq_pop_label
  summarise(count = n(), .groups = "drop")

View(summary_df)

summary_df2 <- dms_maf2$meta$analyses %>% as.data.frame() %>%
  group_by(morphid2, svdq_pop_label) %>%#svdq_pop_label
  summarise(count = n(), .groups = "drop")
View(summary_df2)

# Add a column to indicate the dataset
s3 <- s3 %>% mutate(dataset = "s3")
s4 <- s4 %>% mutate(dataset = "s4")

# Select relevant columns, including `n`
s3_selected <- s3 %>% select(genetic_group, site, Ho, uHe, uFis, loci, dataset)
s4_selected <- s4 %>% select(genetic_group, site, Ho, uHe, uFis, loci, dataset)

# Combine both datasets
combined <- bind_rows(s3_selected, s4_selected)

# Reshape to wide format
wide_result <- combined %>%
  pivot_wider(
    id_cols = c(genetic_group, site),
    names_from = dataset,
    values_from = c(Ho, uHe, uFis, loci)
  )

wide_result <- cbind(wide_result, n=s3$n)
wide_result
write.xlsx(wide_result, file = "LantCama/outputs/LantCama_cluster_stats_80miss_maf0_2.xlsx")

#######

# Add a column to indicate the dataset
s5 <- s5 %>% mutate(dataset = "s5")
s6 <- s6 %>% mutate(dataset = "s6")

# Select relevant columns, including `n`
s2_selected <- s5 %>% select(genetic_group, site, Ho, uHe, uFis, loci, dataset)
s6_selected <- s6 %>% select(genetic_group, site, Ho, uHe, uFis, loci, dataset)

# Combine both datasets
combined2 <- bind_rows(s2_selected, s6_selected)

# Reshape to wide format
wide_result2 <- combined2 %>%
  pivot_wider(
    id_cols = c(genetic_group, site),
    names_from = dataset,
    values_from = c(Ho, uHe, uFis, loci)
  )

wide_result2 <- cbind(wide_result2, n=s5$n)

wide_result3 <- merge(wide_result2, summary_df, by.x='site', by.y='pop_morphEA')

# View the wide-format dataframe
print(wide_result3)

write.xlsx(wide_result3, file = "LantCama/outputs/LantCama_pop_morphEA_stats_80miss_maf0_2.xlsx")

######

fis_maf0 <- individual_F(dms$gt, genetic_group_variable = rep('lantana camara', length(dms$sample_names)), maf = 0, max_missingness = 1)

ho_maf0 <- individual_Ho(dms$gt, genetic_group_variable = rep('lantana camara', length(dms$sample_names)), maf = 0.0, max_missingness = 1)
ho_maf2 <- individual_Ho(dms_maf2$gt, genetic_group_variable = rep('lantana camara', length(dms_maf2$sample_names)), maf = 0, max_missingness = 1)

ho_maf0_df <- data.frame(ho0=as.numeric(ho_maf0), ho2=as.numeric(ho_maf2),dms$meta$analyses) 

# Perorm pairwise comparisons
cm <- compare_means(ho0 ~ morphid2,  data = ho_maf0_df, method = "wilcox.test",p.adjust.method = "fdr") %>% as.data.frame()

cm <- cm[order(cm$group1), ]
my_comparisons <- list()
i=1
signif <- c()
for(row in 1:nrow(cm)){
  if(isTRUE(cm[row, "p.adj"] <= 0.001)){
    my_comparisons[[i]] <- cm[row,c('group1', 'group2')] %>% unlist() %>% as.vector()
    i <- i+1
    # signif <- c(signif, cm[row,'p.adj'] %>% unlist())
    signif <- c(signif, cm[row,'p.signif'] %>% unlist())
  }
}

# start_y_pos <- max(wide_result3$Ho_s2)+0.002
y_positions <- seq(from = max(as.numeric(ho_maf0_df$ho0))+0.001, by = 0.001, length.out = i)
# Create the plot
individual_ho_plot <- ggplot(data = ho_maf0_df, mapping = aes(x = morphid2, y = ho0, fill=morphid2)) +
  geom_boxplot( outliers = TRUE, show.legend=FALSE) +
  scale_fill_manual(values=morphid_colours)+
  # geom_jitter() +
  theme_few() +
  geom_signif(comparisons = my_comparisons, annotations=signif, y_position = y_positions,
              textsize = 3,
              margin_top = 0.001,
              tip_length = 0.001,
              size=0.3, vjust = 0.9,
              color="black")+
  labs(x="Morphotype", y="Individual Ho")
individual_ho_plot


# Perorm pairwise comparisons
cm2 <- compare_means(ho2 ~ morphid2,  data = ho_maf0_df, method = "wilcox.test",p.adjust.method = "fdr") %>% as.data.frame()

cm2 <- cm2[order(cm2$group1), ]
my_comparisons2 <- list()
i2=1
signif2 <- c()
for(row in 1:nrow(cm2)){
  if(isTRUE(cm2[row, "p.adj"] <= 0.001)){
    my_comparisons2[[i2]] <- cm2[row,c('group1', 'group2')] %>% unlist() %>% as.vector()
    i2 <- i2+1
    # signif <- c(signif, cm[row,'p.adj'] %>% unlist())
    signif2 <- c(signif2, cm2[row,'p.signif'] %>% unlist())
  }
}

# start_y_pos <- max(wide_result3$Ho_s2)+0.002
y_positions2 <- seq(from = max(as.numeric(ho_maf0_df$ho2))+0.003, by = 0.003, length.out = i)
# Create the plot
individual_ho2_plot <- ggplot(data = ho_maf0_df, mapping = aes(x = morphid2, y = ho2, fill=morphid2)) +
  geom_boxplot( outliers = TRUE, show.legend=FALSE) +
  scale_fill_manual(values=morphid_colours)+
  # geom_jitter() +
  theme_few() +
  geom_signif(comparisons = my_comparisons2, annotations=signif2, y_position = y_positions2,
              textsize = 3,
              margin_top = 0.001,
              tip_length = 0.001,
              size=0.3, vjust = 0.9,
              color="black")+
  labs(x="Morphotype", y="Individual Ho (MAF 2%)")
individual_ho2_plot

individual_Ho_plots <- ggarrange(individual_ho_plot,individual_ho2_plot, ncol=2, labels="AUTO")

ggsave("LantCama/outputs/LantCama_individual_ho_plots.pdf",
       individual_Ho_plots, width = 22, height = 10, units = "cm", dpi = 600)

ggsave("LantCama/outputs/LantCama_individual_ho_plots.png",
       individual_Ho_plots, width = 22, height = 10, units = "cm", dpi = 600)


####


#### FST ####
# remove sites where n=4
sppop_freq <- as.data.frame(table(dms$meta$site))
not_n1_sites <- as.vector(sppop_freq[sppop_freq$Freq<5,1]) #remove groups where n<=1
not_n1_samples <- dms$sample_names[which(!(dms$meta$site %in% not_n1_sites)& !is.na(dms$meta$site))]
fst_dms <- remove.by.list(dms, not_n1_samples)

length(fst_dms$sample_names)
length(fst_dms$locus_names)

gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$site, RandRbase,species,dataset, maf_val=0.02, miss_val=0.8) #calculates genetic distance
pS        <- population.pw.spatial.dist(fst_dms, fst_dms$meta$site) #calculates geographic distance between populations

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

# adding metadata for sites
Fst_sig2 <- merge(Fst_sig, distinct(m2[,c("site","morphid2")]), by.x="Var1", by.y="site", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(m2[,c("site","morphid2")]), by.x="Var2", by.y="site", all.y=FALSE)
Fst_sig2$same_morphid2 <- ifelse(Fst_sig2$morphid2.x == Fst_sig2$morphid2.y, "Intra-morph", "Inter-morph")

library(ggforce)
fstp1 <- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_morphid2))+geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  # facet_zoom(x=Geo_dist2<25, zoom.size=1)+
  theme_bw()+
  geom_hline(yintercept = 0.3, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1

# 
# 
# # adding metadata for sites
# Fst_sig2 <- merge(Fst_sig, distinct(m2[,c("site","svdq_pop_label")]), by.x="Var1", by.y="site", all.y=FALSE)
# Fst_sig2 <- merge(Fst_sig2, distinct(m2[,c("site","svdq_pop_label")]), by.x="Var2", by.y="site", all.y=FALSE)
# Fst_sig2$same_svdq_pop_label <- ifelse(Fst_sig2$svdq_pop_label.x == Fst_sig2$svdq_pop_label.y, "Intra-morph", "Inter-morph")
# 
# library(ggforce)
# fstp1 <- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_svdq_pop_label))+geom_point(size=1, alpha=0.3)+
#   labs(x="Distance (km)", y="FST", colour="Comparison")+
#   geom_hline(yintercept = 0.3, linetype="dotted")+
#   # facet_zoom(x=Geo_dist2<25, zoom.size=1)+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
# fstp1

# ggsave("BossFrag/outputs/paper/supfig2_BossFrag_manning_fst.pdf",
#        fstp1, width = 15, height = 15, units = "cm", dpi=600)

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

agg <- unique(m2[, c("site", "morphid2",'national2')]) # create aggregated df of pop_largeecies and site
mat2 <- merge(mat2, agg, by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

row_group_ann <- rowAnnotation(Morphotype = mat2$morphid2,
                               col=list(Morphotype=morphid_colours),
                               na_col="white",
                               annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                              title_gp=gpar(fontsize=10)),
                               annotation_name_gp = gpar(fontsize = 0),
                               annotation_name_side="top")



bottom_group_ann <- HeatmapAnnotation(Morphotype = mat2$morphid2, col = list(Morphotype = morphid_colours),
                                      annotation_name_gp = gpar(fontsize = 0),
                                      annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
                                                                     title_gp=gpar(fontsize=10)),
                                      annotation_name_side="left",
                                      na_col = "white")



row_group_ann2 <- rowAnnotation(Country = mat2$national2,
                               col=list(Country=nation_colours),
                               na_col="white",
                               annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                              title_gp=gpar(fontsize=10)),
                               annotation_name_gp = gpar(fontsize = 0),
                               annotation_name_side="top")



bottom_group_ann2 <- HeatmapAnnotation(Country = mat2$national2, col = list(Country = nation_colours),
                                      annotation_name_gp = gpar(fontsize = 0),
                                      annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
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
                })

gene_width <- nrow(mat2)*unit(6, "mm")

# draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

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


#### pca ####

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
  guides(colour = guide_legend(title.position = "top"))+#+
  scale_colour_manual(values=morphid_colours)
# scale_shape_manual(values=cluster2ecies_shapes)
pca_plot1

pca_plot2 <- ggplot(g_pca_df2, aes(x=PC3, y=PC4, colour=morphid2, shape=national2))+ xlab(pcnames[3])+ylab(pcnames[4])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)
# scale_colour_manual(values=cluster2ecies_colours)+
# scale_shape_manual(values=cluster2ecies_shapes)

pca_plot3 <- ggplot(g_pca_df2, aes(x=PC5, y=PC6, colour=morphid2, shape=national2))+ xlab(pcnames[5])+ylab(pcnames[6])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)
# scale_colour_manual(values=cluster2ecies_colours)+
# scale_shape_manual(values=cluster2ecies_shapes)

all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3, labels=c("A","B","C"),
                            common.legend = TRUE, ncol=3, legend = "bottom")
all3_pca_plots
