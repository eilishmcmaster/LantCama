library(dplyr)
library(ggpubr)
library(popkin)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ggthemes)
library(ggplot2)
library(ggpubr)
### Import data ###

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/kin_compare_functions.R?raw=TRUE")


# meta <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/Lcam_DLan23-8067_meta.RRv0002.csv", sep=",")
meta <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/meta_targetid.csv")

#### Make colour palette ####
cluster_colours <- scales::hue_pal()(length(unique(meta$cluster)))
names(cluster_colours) <- unique(meta$cluster)[order(unique(meta$cluster))]

cluster_shapes <- 1:(length(unique(meta$cluster)))
names(cluster_shapes) <- unique(meta$cluster)[order(unique(meta$cluster))]


#### PCA of all lineages ####
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')

genotype_matrix <- vcfR::vcfR2genind(vcfraw, return.alleles = TRUE)

genotype_matrix$tab[1:10,1:10]

gt <- genotype_matrix$tab

loc_missing <- colMeans(is.na(gt))
hist(loc_missing)
gt <- gt[,loc_missing<=0.2]
ncol(gt)

sample_missing <- rowMeans(is.na(gt))
hist(sample_missing)
gt <- gt[sample_missing<=0.2,]
nrow(gt)
#### need to remove high missing samples



gen_d5 <- new("genlight", gt) #convert df to genlight object for glPca function
unique(ploidy(gen_d5)) # adegenet knows that its tetraploid 

gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] #extract PCs 
g_pca_df2 <- merge(g_pca_df, meta, by.x=0, by.y="targetid", all.y=FALSE, all.x=FALSE) # add metadata 


pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)")


pca_plot <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=cluster, shape=cluster))+ 
  geom_point(alpha=0.7)+theme_few()+xlab(pcnames[1])+ylab(pcnames[2])+
  labs(colour="Cluster", shape="Cluster")+
  scale_colour_manual(values=cluster_colours)+
  scale_shape_manual(values=cluster_shapes)+  guides(
    shape = guide_legend(keywidth = 1, keyheight = 1),
    color = guide_legend(keywidth = 1, keyheight = 1)
  )


pca_plot

##### PCA plot of EACP only ####

vcfraw_eacp <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered_eacp.vcf')

genotype_matrix_eacp <- vcfR::vcfR2genind(vcfraw_eacp, return.alleles = FALSE)

genotype_matrix_eacp$tab[1:10, 1:10]

gt_eacp <- genotype_matrix_eacp$tab

loc_missing_eacp <- colMeans(is.na(gt_eacp))
hist(loc_missing_eacp)
gt_eacp <- gt_eacp[, loc_missing_eacp <= 0.2]
ncol(gt_eacp)
sample_missing_eacp <- rowMeans(is.na(gt_eacp))
hist(sample_missing_eacp)
gt_eacp <- gt_eacp[sample_missing_eacp<=0.2,]
nrow(gt_eacp)

gen_d5_eacp <- new("genlight", gt_eacp) #convert df to genlight object for glPca function
unique(ploidy(gen_d5_eacp)) # adegenet knows that its tetraploid 

gen_pca_eacp <- glPca(gen_d5_eacp, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df_eacp <- gen_pca_eacp[["scores"]] #extract PCs 
g_pca_df2_eacp <- merge(g_pca_df_eacp, meta, by.x = 0, by.y = "targetid", all.y = FALSE, all.x = FALSE) # add metadata 

pcnames_eacp <- paste0(colnames(g_pca_df_eacp), " (",
                       paste(round(gen_pca_eacp[["eig"]][1:5] / sum(gen_pca_eacp[["eig"]]) * 100, 2)),
                       "%)")

pca_plot_eacp <- ggplot(g_pca_df2_eacp, aes(x = PC1, y = PC2, colour = cluster, shape=cluster)) + 
  geom_point() + theme_few() + xlab(pcnames_eacp[1]) + ylab(pcnames_eacp[2]) +
  labs(colour="Cluster", shape="Cluster")+
  scale_colour_manual(values=cluster_colours)+
  scale_shape_manual(values=cluster_shapes)+
  theme(legend.spacing.x = unit(0, 'cm'))

pca_plot_eacp
pca_plots_comb <- ggarrange(pca_plot, pca_plot_eacp, common.legend = TRUE, labels=c("A","B"),legend="right")

ggsave("LantCama/outputs/plots/vcf_plots/vcf_pca_plots.png", 
       plot = pca_plots_comb, width = 180, height = 80, dpi = 300, units = "mm")


#### Ploidy investigation by readcount ####
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')
# vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/20230920_lantana_dp_qual.vcf')
# dp <- extract.gt(vcfraw, element = 'DP') #allele depth

ad <- extract.gt(vcfraw, element = 'AD') #allele depth 
# gta <- extract.gt(vcfraw, element = 'GT')
# hets <- is_het(gta)
#remove homos
# is.na( ad[ !hets ] ) <- TRUE

c1o <- masplit(ad, record = 1)
c1 <- c1o
c1 <- c1[,colnames(c1) %in% rownames(gt)]
matching_targetids <- meta$targetid[meta$targetid %in% colnames(c1)]
c1 <- c1[, colnames(c1) %in% matching_targetids]
matching_sample <- meta$sample[meta$targetid %in% matching_targetids]
colnames(c1) <- matching_sample

c2o <- masplit(ad, record = 2)
c2 <- c2o
c2 <- c2[,colnames(c2) %in% rownames(gt)]
matching_targetids <- meta$targetid[meta$targetid %in% colnames(c2)]
c2 <- c2[, colnames(c2) %in% matching_targetids]
matching_sample <- meta$sample[meta$targetid %in% matching_targetids]
colnames(c2) <- matching_sample

vcf_counts <- list(c1=c1, c2=c2)


vcf_test <- read_histogram_function2(meta=meta, vcf_counts,
                                     run_quantile = TRUE, min_quantile=0.1, max_quantile=0.9,
                                 min_depth=10, species_col="cluster")





#### whole species histograms ####

vcf_z <- whole_sp_plots(vcf_test,  c("eacp", "eawt", "per1"), NULL)
vcf_sp_hist_plots <- ggarrange(vcf_z[[1]],vcf_z[[2]],vcf_z[[3]], align="hv", ncol=3,
                           labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
vcf_sp_hist_plots

ggsave("LantCama/outputs/plots/vcf_plots/vcf_species_ploidy_hist.png", plot = vcf_sp_hist_plots, width = 150, height = 60, dpi = 300, units = "mm")

##### save individual plots to pdf to find best examples #####
specific_sample_plots2 <- function(data){
  plots <- list()
  samples <- rownames(data)
  for(i in samples){
    print(i)
    x <- data[i,] %>% as.numeric()
    z <- data.frame(x=x)
    
    p <- ggplot(z, aes(x)) +
      geom_histogram(bins = 30, fill = "lightblue", color="black", size=0.2) + #usually use 50 bins
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4 ", "1/3", "1/2", "2/3", " 3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(i))+
      theme_few()+
      # scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black"), angle=90),
            axis.text.y = element_text(size=8),
            title = element_text(size=8),
            plot.margin = margin(0,0,0,0))
    
    
    plots[[paste0(i)]] <- p
    rm(x)
  }
  return(plots)
}


eacp_all_plots <- specific_sample_plots2(vcf_test$eacp)

pdf("/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/plots/vcf_plots/eacp_allplots.pdf",width = 3, height = 3, onefile = TRUE)
for(i in 1:length(eacp_all_plots)){
  tplot <- eacp_all_plots[[i]]
  print(tplot)
}
dev.off()

eacp_good_samples <- c("NSW1078961", "NSW1152126","NSW1089413","NSW1089204","NSW1095153","NSW1096750.1")

per1_all_plots <- specific_sample_plots2(vcf_test$per1)
pdf("/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/plots/vcf_plots/per1_allplots.pdf",width = 3, height = 3, onefile = TRUE)
for(i in 1:length(per1_all_plots)){
  tplot <- per1_all_plots[[i]]
  print(tplot)
}
dev.off()

per1_good_samples <- c("NSW1089128","NSW1150350","NSW1150436","NSW1158956","NSW1161295")


eawt_all_plots <- specific_sample_plots2(vcf_test$eawt)
pdf("/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/plots/vcf_plots/eawt_allplots.pdf",width = 3, height = 3, onefile = TRUE)
for(i in 1:length(eawt_all_plots)){
  tplot <- eawt_all_plots[[i]]
  print(tplot)
}
dev.off()

eawt_good_samples <- c("NSW1084601","NSW1096829","NSW1089497")

### individual esample histograms ####

vcf_eacp_samples <- specific_sample_plots(vcf_test$eacp,
                                          c("NSW1095153", "NSW1152126","NSW1089413"))

vcf_eawt_samples <- specific_sample_plots(vcf_test$eawt,
                                          c("NSW1084601","NSW1096829","NSW1089497"))

vcf_per1_samples <- specific_sample_plots(vcf_test$per1,
                                      c("NSW1089128","NSW1150350","NSW1150436"))

vcf_all_hist <- ggarrange(vcf_eacp_samples[[1]],vcf_eacp_samples[[2]],vcf_eacp_samples[[3]],
                      vcf_eawt_samples[[1]],vcf_eawt_samples[[2]],vcf_eawt_samples[[3]],
                      vcf_per1_samples[[1]],vcf_per1_samples[[2]],vcf_per1_samples[[3]],
                      align="hv", ncol=3, nrow=3,
                      labels=c("A","","","B","","","C"),
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL))%>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count"
  )

# vcf_all_hist

ggsave("LantCama/outputs/plots/vcf_plots/vcf_all_ploidy_hist2.png", plot = vcf_all_hist, width = 190, height = 170, dpi = 300, units = "mm")

###
# https://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf

ad <- extract.gt(vcfraw, element = 'AD') #allele depth 
gt <- extract.gt(vcfraw, element = 'GT')
hets <- is_het(gt)

#take a look
gt[1:10,1:10]

#remove homos
is.na( ad[ !hets ] ) <- TRUE

allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)

hist(ad2, breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1, breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)

axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# get 15% and 95% quantiles for most abundant alleles
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)

sums[,1]
# most abundant allele hist with 15% and 95% quantiles for locus 1
tmp <- allele1[,1]
hist(tmp, breaks=seq(0,100,by=1), col="#808080", main = "P17777us22")
abline(v=sums[,1], col=2, lwd=2)
# second most abundant allele hist with 15% and 95% quantiles for locus 1
tmp2 <- allele2[,1]
hist(tmp2, breaks=seq(0,100,by=1), col="#808080", main = "P17777us22")
abline(v=sums[,1], col=2, lwd=2)

# boxplot the allele depth on the first 10 loci 
boxplot(allele1[,1:10], las=3)


# FILTER ALL DATA TO REMOVE RESULTS THAT DONT MEET ALLELE DEPTH REQUIREMENTS
# Assign the value of vcfraw to a new variable called vcf
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')
vcf <- vcfraw
# Extract the genotype information for allele depth (AD) and store it in the ad variable
ad <- extract.gt(vcf, element = 'AD')

# Split allele depth data into allele 1 and allele2
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2) 

# Calculate the 15th and 95th percentiles for each column (variant) in allele1
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)

# Subtract the 15th percentile values from each column in allele1 and store the result in dp2
dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[1,])
# Set genotype values to NA for positions where dp2 is less than 0 and the genotype is not already NA
vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# Subtract the 95th percentile values from each column in allele1 and store the result in dp2
dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[2,])
# Set genotype values to NA for positions where dp2 is greater than 0
vcf@gt[,-1][dp2 > 0] <- NA
# Subtract the 15th percentile values from each column in allele2 and store the result in dp2
dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[1,])
# Set genotype values to NA for positions where dp2 is less than 0 and the genotype is not already NA
vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# Subtract the 95th percentile values from each column in allele2 and store the result in dp2
dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[2,])
# Set genotype values to NA for positions where dp2 is greater than 0
vcf@gt[,-1][dp2 > 0] <- NA


ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
boxplot(allele1[,1:10], las=3)


gt <- extract.gt(vcf, element = 'GT')
hets <- is_het(gt)
is.na( ad[ !hets ] ) <- TRUE

allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

# Set the file name for the PNG
png_file <- "/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/plots/vcf_plots/all_sp_filtered_knausb_method_plot.png"
# Open a PNG graphics device
png(png_file, width = 800, height = 600)  # Adjust width and height as needed
# Create the histogram plot for ad2
hist(ad2[,1], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
# Add the histogram plot for ad1
hist(ad1[,1], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
# Customize the x-axis labels
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
# Close the graphics device and save the plot to the PNG file
dev.off()


# Eilish plot method
# vcf_all_counts <- list()
# vcf_all_counts$c1 <- allele1
# vcf_all_counts$c2 <- allele2
# # need to get the meta meta with targetids
# test <- read_histogram_function2(meta=meta, counts=counts2,
#                                  min_depth=10, min_quantile=0.05, max_quantile=0.95, species_col="sp")
###



###############
