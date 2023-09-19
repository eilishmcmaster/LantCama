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


m2 <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/Lcam_DLan23-8067_meta.RRv0002.csv", sep=",")
c2 <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/meta_targetid.csv")

#### Make colour palette ####
cluster_colours <- scales::hue_pal()(length(unique(c2$cluster)))
names(cluster_colours) <- unique(c2$cluster)[order(unique(c2$cluster))]

cluster_shapes <- 1:(length(unique(c2$cluster)))
names(cluster_shapes) <- unique(c2$cluster)[order(unique(c2$cluster))]


#### PCA of all lineages ####
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')

genotype_matrix <- vcfR::vcfR2genind(vcfraw, return.alleles = FALSE)

genotype_matrix$tab[1:10,1:10]

gt <- genotype_matrix$tab

loc_missing <- colMeans(is.na(gt))
hist(loc_missing)
gt <- gt[,loc_missing<=0.3]
ncol(gt)
sample_missing <- rowMeans(is.na(gt))
hist(sample_missing)
gt <- gt[sample_missing<=0.3,]
nrow(gt)
#### need to remove high missing samples

gen_d5 <- new("genlight", gt) #convert df to genlight object for glPca function
unique(ploidy(gen_d5)) # adegenet knows that its tetraploid 

gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] #extract PCs 
g_pca_df2 <- merge(g_pca_df, c2, by.x=0, by.y="targetid", all.y=FALSE, all.x=FALSE) # add metadata 


pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)")


pca_plot <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=cluster, shape=cluster))+ 
  geom_point(alpha=0.7)+theme_few()+xlab(pcnames[1])+ylab(pcnames[2])+
  labs(colour="Genetic group", shape="Genetic group")+
  scale_colour_manual(values=cluster_colours)+
  scale_shape_manual(values=cluster_shapes)

pca_plot

##### PCA plot of EACP only ####

vcfraw_eacp <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered_eacp.vcf')

genotype_matrix_eacp <- vcfR::vcfR2genind(vcfraw_eacp, return.alleles = FALSE)

genotype_matrix_eacp$tab[1:10, 1:10]

gt_eacp <- genotype_matrix_eacp$tab

loc_missing_eacp <- colMeans(is.na(gt_eacp))
hist(loc_missing_eacp)
gt_eacp <- gt_eacp[, loc_missing_eacp <= 0.3]
ncol(gt_eacp)
sample_missing_eacp <- rowMeans(is.na(gt_eacp))
hist(sample_missing_eacp)
gt_eacp <- gt_eacp[sample_missing_eacp<=0.3,]
nrow(gt_eacp)

gen_d5_eacp <- new("genlight", gt_eacp) #convert df to genlight object for glPca function
unique(ploidy(gen_d5_eacp)) # adegenet knows that its tetraploid 

gen_pca_eacp <- glPca(gen_d5_eacp, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df_eacp <- gen_pca_eacp[["scores"]] #extract PCs 
g_pca_df2_eacp <- merge(g_pca_df_eacp, c2, by.x = 0, by.y = "targetid", all.y = FALSE, all.x = FALSE) # add metadata 

pcnames_eacp <- paste0(colnames(g_pca_df_eacp), " (",
                       paste(round(gen_pca_eacp[["eig"]][1:5] / sum(gen_pca_eacp[["eig"]]) * 100, 2)),
                       "%)")

pca_plot_eacp <- ggplot(g_pca_df2_eacp, aes(x = PC1, y = PC2, colour = cluster, shape=cluster)) + 
  geom_point() + theme_few() + xlab(pcnames_eacp[1]) + ylab(pcnames_eacp[2]) +
  labs(colour="Genetic group", shape="Genetic group")+
  scale_colour_manual(values=cluster_colours)+
  scale_shape_manual(values=cluster_shapes)

pca_plot_eacp
ggarrange(pca_plot, pca_plot_eacp, common.legend = TRUE, legend="right")


#### Ploidy investigation by readcount ####
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
# # need to get the m2 meta with targetids
# test <- read_histogram_function2(meta=m2, counts=counts2,
#                                  min_depth=10, min_quantile=0.05, max_quantile=0.95, species_col="sp")
###

### FILTER LOCI BY MISSINGNESS ####
gt <- extract.gt(vcf, element = 'GT')

gt1 <- masplit(gt, record = 1, delim="/")
gt2 <- masplit(gt, record = 2, delim="/")
gt3 <- masplit(gt, record = 3, delim="/")
gt4 <- masplit(gt, record = 4, delim="/")

gta <- (gt1 + gt2 + gt3 + gt4)

loci_missingness <- rowMeans(is.na(gta))
hist(loci_missingness)
gta2 <- gta[loci_missingness<=0.3,]
nrow(gta2)
loci_missingness2 <- rowMeans(is.na(gta2))
hist(loci_missingness2)
hist(colMeans(is.na(gta2)))
