
library(openxlsx) #reading and writing xlsx
library(stringr) #replacing strings while data wrangling 
library(dplyr) #data wrangling 
library(data.table) #data wrangling 
library(ggfortify) # calculating glm and pca
library(ggpubr) # ggarrange
library(pracma) # maths
library(RRtools) #Jason's package for dart data
library(ggthemes) #themes for ggplots
library(RColorBrewer) #used for making colour scemes for plots
library(ozmaps) #draws australia coastlines and state boundaries
library(adegenet) #essential for processing dart data
library(ggrepel) #used for plotting labels on ggplote
library(openxlsx)
library(readxl)


topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory 
species <- "LantCama" #species name
dataset <- "DLan22-7500" #dart order
missingness <- 0.3

# source my custom functions
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")

setwd("/Users/eilishmcmaster/Documents/LantCama")

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")



#######################
# this filters the readcounts by minimum and maximum quantiles which is done in th vcf method as well
filter_count_by_quantiles <- function(count_df, min_quantile, max_quantile){
  sums <- apply(count_df, MARGIN=2, quantile, probs=c(min_quantile, max_quantile), na.rm=TRUE)
  dp1 <- sweep(count_df, MARGIN=2, FUN = "-", sums[1,])
  dp2 <- sweep(count_df, MARGIN=2, FUN = "-", sums[2,])
  out_df <- count_df
  out_df[dp1<0 | dp2>0] <- 0
  return(out_df)
}



read_histogram_function2 <- function(meta, counts, min_depth,run_quantile=NULL, min_quantile=NULL, max_quantile=NULL, species_col, dms=NULL, remove_by_dms=NULL) {
  c1 <- counts$c1
  c2 <- counts$c2
  
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
  # filter the reads
  # filter by quantiles
  if(isTRUE(run_quantile)){
    sums <- apply(c1, MARGIN=2, quantile, probs=c(min_quantile, max_quantile), na.rm=TRUE)
    
    dp1 <- sweep(c1, MARGIN=2, FUN = "-", sums[1,])
    dp2 <- sweep(c1, MARGIN=2, FUN = "-", sums[2,])
    c1[dp1<0 | dp2>0] <- NA
    
    dp12 <- sweep(c2, MARGIN=2, FUN = "-", sums[1,])
    dp22 <- sweep(c2, MARGIN=2, FUN = "-", sums[2,])
    c2[dp12<0 | dp22>0] <- NA
  }
  
  c1[is.na(c1)] <- 0
  c2[is.na(c2)] <- 0
  

  #filter by total number of reads (sometimes lower quantile is 0)
  combined_reads <- c1 + c2
  c1[combined_reads < min_depth] <- 0
  c2[combined_reads < min_depth] <- 0
  combined_reads <- c1 + c2
  
  
  # get the proportions for all (rows are samples)
  c3_min <- pmin(t(c1), t(c2), na.rm = TRUE) / t(combined_reads)
  c3_min[is.infinite(c3_min)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_min[is.nan(c3_min)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  c3_max <- pmax(t(c1), t(c2), na.rm = TRUE) / t(combined_reads)
  c3_max[is.infinite(c3_max)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_max[is.nan(c3_max)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  # removes full homozygotes otherwise theres waaaay too many to plot
  c3_min[c3_min==0] <- NA
  c3_min[c3_min==1] <- NA
  c3_max[c3_max==0] <- NA
  c3_max[c3_max==1] <- NA
  
  c3 <- cbind(c3_min, c3_max) 
  
  out_data  <- list() # place to put the data 
  
  for (i in seq_along(species)) {
    print(paste("Running", species[i], "now"))
    samples <- meta$sample[meta[[species_col]] == species[i]] # get the NSW ID for that species samples
    if(isTRUE(remove_by_dms)){
      samples <- samples[samples %in% dms$sample_names]
    }
    
    c3_species <- c3[row.names(c3) %in% samples, ] # get the readcount df with those samples
    
    if (class(c3_species) %in% "array" || class(c3_species) %in% "matrix") {
      c3_species <- c3_species[, colSums(!is.na(c3_species) & c3_species != "") > 0] 
      out_data[[paste0(species[i])]] <- data.frame(c3_species)
    } else {
      c3_species <- c3_species[!is.na(c3_species)] # remove empty columns
      out_data[[paste0(species[i])]] <- as.vector(c3_species)
    }
  }
  
  return(out_data)
}


whole_sp_plots <- function(data, species, max){
  plots <- list()
  for(i in seq_along(species)){
    x <- na.omit(unlist(list(data[[species[i]]])))
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 30, fill = "lightpink", color="black", size=0.2) + 
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4", "1/3", "1/2", "2/3", "3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(species[i]))+
      theme_few()+
      # scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black"), angle=90),
            axis.text.y = element_text(size=8),
            title = element_text(size=8, face="italic"),
            plot.margin = margin(0, 0, 0, 0))
    
    plots[[i]] <- p
  }
  return(plots)
}

specific_sample_plots <- function(data, samples){
  plots <- list()
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

#### DArT ####
###### Ploidy investigation by DArT ####
counts2 <- read_dart_counts_csv_faster('LantCama/dart_raw/Report_DLan22-7500_3_moreOrders_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=0, 
                                       minGenotypeCount=0)

m2 <- custom.read(species, dataset) #read custom metadata csv

test <- read_histogram_function2(meta=m2, counts=counts2,run_quantile = FALSE,min_quantile=0.1, max_quantile=0.9,
                                 min_depth=10,  species_col="sp") #min_quantile=0.05, max_quantile=0.95,

####### Whole cluster histograms #####
z <- whole_sp_plots(test,  c("eacp", "eawt", "per1"), NULL)
sp_hist_plots <- ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3,
                           labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
sp_hist_plots

ggsave("LantCama/outputs/plots/dart_species_ploidy_hist.png", plot = sp_hist_plots, width = 150, height = 60, dpi = 300, units = "mm")


####### Specific example histograms #####
eacp_samples <- specific_sample_plots(test$eacp,
                                          c("NSW1095153", "NSW1152126","NSW1089413"))

eawt_samples <- specific_sample_plots(test$eawt,
                                          c("NSW1084601","NSW1096829","NSW1089497"))

per1_samples <- specific_sample_plots(test$per1,
                                          c("NSW1089128","NSW1150350","NSW1150436"))

all_hist <- ggarrange(eacp_samples[[1]],eacp_samples[[2]],eacp_samples[[3]],
                      eawt_samples[[1]],eawt_samples[[2]],eawt_samples[[3]],
                      per1_samples[[1]],per1_samples[[2]],per1_samples[[3]],
                      align="hv", ncol=3, nrow=3,
                      labels=c("A","","","B","","","C"),
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL))%>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count"
  )

# all_hist

ggsave("LantCama/outputs/plots/dart_all_ploidy_hist.png", plot = all_hist, width = 190, height = 170, dpi = 300, units = "mm")

#### VCF #### 

##### Setup ####
# meta <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/Lcam_DLan23-8067_meta.RRv0002.csv", sep=",")
meta <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/meta_targetid.csv")

###### Make colour palette ####
cluster_colours <- scales::hue_pal()(length(unique(meta$cluster)))
names(cluster_colours) <- unique(meta$cluster)[order(unique(meta$cluster))]

cluster_shapes <- 1:(length(unique(meta$cluster)))
names(cluster_shapes) <- unique(meta$cluster)[order(unique(meta$cluster))]


##### PCA of all lineages ####
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')

# convert to genind
genotype_matrix <- vcfR::vcfR2genind(vcfraw, return.alleles = TRUE)

gt <- genotype_matrix$tab

# remove high missing loci
loc_missing <- colMeans(is.na(gt))
hist(loc_missing)
gt <- gt[,loc_missing<=0.2]
ncol(gt)

# remove high missing samples
sample_missing <- rowMeans(is.na(gt))
hist(sample_missing)
gt <- gt[sample_missing<=0.2,]
nrow(gt)

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




##### Ploidy investigation by VCF ####
vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/lantana_full_filtered.vcf')

# get allele depth
ad <- extract.gt(vcfraw, element = 'AD') #allele depth 

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



######## Whole cluster histograms ####

vcf_z <- whole_sp_plots(vcf_test,  c("eacp", "eawt", "per1"), NULL)
vcf_sp_hist_plots <- ggarrange(vcf_z[[1]],vcf_z[[2]],vcf_z[[3]], align="hv", ncol=3,
                               labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
vcf_sp_hist_plots

ggsave("LantCama/outputs/plots/vcf_plots/vcf_species_ploidy_hist.png", plot = vcf_sp_hist_plots, width = 150, height = 60, dpi = 300, units = "mm")

####### Specific example histograms ####

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


ggsave("LantCama/outputs/plots/vcf_plots/vcf_all_ploidy_hist2.png", plot = vcf_all_hist, width = 190, height = 170, dpi = 300, units = "mm")


