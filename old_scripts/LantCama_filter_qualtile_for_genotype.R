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


genotype_matrix <- vcfR::vcfR2genind(vcf, return.alleles = TRUE)


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

