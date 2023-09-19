
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

counts2 <- read_dart_counts_csv_faster('LantCama/dart_raw/Report_DLan22-7500_3_moreOrders_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=0, 
                                       minGenotypeCount=0)

library(readxl)
m2 <- custom.read(species, dataset) #read custom metadata csv

#######################
# this filters the readcounts by minimum and maximum quantiles which is done in th vcf method as well
filter_count_by_quantiles <- function(count_df, min_quantile, max_quantile){
  quantiles <- apply(count_df, MARGIN=2, quantile, probs=c(min_quantile, max_quantile), na.rm=TRUE)
  dp1 <- sweep(count_df, MARGIN=2, FUN = "-", sums[1,])
  dp2 <- sweep(count_df, MARGIN=2, FUN = "-", sums[2,])
  out_df <- count_df
  out_df[dp1<0 | dp2>0] <- 0
  return(out_df)
}


read_histogram_function2 <- function(meta, counts, min_depth, min_quantile, max_quantile, species_col, dms=NULL, remove_by_dms=NULL) {
  c1 <- counts$c1
  c2 <- counts$c2
  
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
  # anything that is NA must be 0 for the dividing etc
  c1[is.na(c1)] <- 0
  c2[is.na(c2)] <- 0
  
  # filter the reads
  # filter by quantiles
  c1 <- filter_count_by_quantiles(c1, min_quantile, max_quantile)
  c2 <- filter_count_by_quantiles(c2, min_quantile, max_quantile)
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


test <- read_histogram_function2(meta=m2, counts=counts2,
                                 min_depth=10, min_quantile=0.15, max_quantile=0.95, species_col="sp")


####

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

# counts3 <- count_subsetter(dms, counts2, 0)


z <- whole_sp_plots(test,  c("eacp", "eawt", "per1"), NULL)
sp_hist_plots <- ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3,
                           labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
sp_hist_plots

z[[1]]
# sp_hist_plots
ggsave("LantCama/outputs/plots/species_ploidy_hist2.png", plot = sp_hist_plots, width = 150, height = 60, dpi = 300, units = "mm")


###
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
##
eacp_samples <- specific_sample_plots(test$eacp,
                                      c("NSW1089413","NSW1096776","NSW1095152"))#c(160,160,160)) # for 50 breaks

eawt_samples <- specific_sample_plots(test$eawt,
                                      c("NSW1084671","NSW1084666","NSW1095126"))#c(250,250,250))

per1_samples <- specific_sample_plots(test$per1,
                                      c("NSW1158953","NSW1150367","NSW1161296"))#c(50,50,50))

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

all_hist

ggsave("LantCama/outputs/plots/all_ploidy_hist.png", plot = all_hist, width = 190, height = 170, dpi = 300, units = "mm")


#####################

