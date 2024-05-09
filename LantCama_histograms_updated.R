library(patchwork)
species <- "LantCama" #species name
dataset <- "DLan23-8067" #dart order


#### DArT ####
###### Ploidy investigation by DArT ####
counts2 <- read_dart_counts_csv_faster('LantCama/dart_raw/Report_DLan22-7500_3_moreOrders_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=0, 
                                       minGenotypeCount=0)

m2 <- custom.read(species, dataset) %>% as.data.frame()#read custom metadata csv

test <- read_histogram_function2(meta=m2, counts=counts2,run_quantile = FALSE,min_quantile=0.1, max_quantile=0.9,
                                 min_depth=10,  species_col="cluster_histogram") #min_quantile=0.05, max_quantile=0.95,

####### Specific example histograms #####
test2 <- read_histogram_function2(meta=m2, counts=counts2,run_quantile = FALSE,min_quantile=0.1, max_quantile=0.9,
                                 min_depth=10,  species_col="site_histogram") #min_quantile=0.05, max_quantile=0.95,

z <- whole_sp_plots(test2,  names(test2), NULL)
sp_hist_plots <- ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3,
                           labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
sp_hist_plots

ggsave("LantCama/outputs/plots/dart_species_ploidy_hist.png", plot = sp_hist_plots, width = 250, height = 100, dpi = 300, units = "mm")


coot_samples <- specific_sample_plots(test2$`Mt Coot-tha`,rownames(test2$`Mt Coot-tha`))
coot_plots <- wrap_plots(coot_samples, nrow = 2, ncol = 3)

wilb_samples <- specific_sample_plots(test2$Wilberforce,rownames(test2$Wilberforce))
wilb_plots <- wrap_plots(wilb_samples, nrow = 2, ncol = 3)

norah_samples <- specific_sample_plots(test2$`Norah Head`,rownames(test2$`Norah Head`))
norah_plots <- wrap_plots(norah_samples, nrow = 2, ncol = 3)

all_hist <- ggarrange(coot_plots,wilb_plots,norah_plots,ncol=1, labels="AUTO")

ggsave("LantCama/outputs/plots/dart_all_ploidy_hist.png", plot = all_hist, width = 250, height = 300, dpi = 300, units = "mm")
