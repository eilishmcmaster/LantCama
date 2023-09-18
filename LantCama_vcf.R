library(dplyr)
library(ggpubr)
library(popkin)
library(SNPRelate)
library(vcfR)

### Import data ###

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/kin_compare_functions.R?raw=TRUE")


m2 <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/Lcam_DLan23-8067_meta.RRv0002.csv", sep=",")

# vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/raw_incomp.mac3bcf.vcf')
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered_eacp.vcf')

genotype_matrix <- vcfR::vcfR2genind(vcfraw, return.alleles = FALSE)

loc_missing <- colMeans(is.na(genotype_matrix$tab))

gen_d5 <- new("genlight", genotype_matrix$tab) #convert df to genlight object for glPca function

gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] #extract PCs 
# g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # add metadata 


pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes


pca_plot2 <- ggplot(g_pca_df, aes(x=PC1, y=PC2))+ 
  geom_point()+theme_few()+xlab(pcnames[1])+ylab(pcnames[2])
# scale_colour_manual(values=gro)
pca_plot2

#####


#https://knausb.github.io/vcfR_documentation/determining_ploidy_1.html


vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')

genotype_matrix <- vcfR::vcfR2genind(vcfraw, return.alleles = TRUE)

## each strand is put into a different column????
loc_missing <- colMeans(is.na(genotype_matrix$tab))
hist(loc_missing)
gt <- genotype_matrix$tab[,which(loc_missing<=0.3)]
ncol(gt)

gen_d5 <- new("genlight", gt) #convert df to genlight object for glPca function
ploidy(gen_d5)

gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] #extract PCs 
# g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # add metadata 


pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes


pca_plot2 <- ggplot(g_pca_df, aes(x=PC1, y=PC2))+ 
  geom_point()+theme_few()+xlab(pcnames[1])+ylab(pcnames[2])
# scale_colour_manual(values=gro)
pca_plot2

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


gt_tidy <- extract_gt_tidy(vcfraw)
