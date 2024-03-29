library(readr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
#https://speciationgenomics.github.io/filtering_vcfs/
# LantCama Vcf qc
setwd("/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf_qc/")


custom_theme <- theme(
  axis.text = element_text(size = 8),    # Adjust the size for axis labels
  axis.title = element_text(size = 10)   # Adjust the size for axis titles
)

# variant quality 
var_qual <- read_delim("variant_random_subset.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
# var_qual <- read_delim("lantana_filtered.lqual", delim = "\t",
#                        col_names = c("chr", "pos", "qual"), skip = 1)

var_qc_stats <- summary(var_qual$qual)

var_qc_plot <- ggplot(var_qual, aes(qual)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins=100)+ 
  ylab("count")+
  theme_few()+xlim(1,1000)+geom_vline(xintercept = 20, color="red", linetype="dotted")

var_qc_plot2 <- var_qc_plot+custom_theme+
  annotate("text", size=2, color="red", x=Inf, y=Inf, hjust=1, vjust=1, label = paste("min:", round(var_qc_stats[1],2),
                                                                         "\nQ1:", round(var_qc_stats[2],2),
                                                                         "\nmedian:", round(var_qc_stats[3],2),
                                                                         "\nmean:", round(var_qc_stats[4],2),
                                                                         "\nQ3:", round(var_qc_stats[5],2),
                                                                         "\nmax:", round(var_qc_stats[6],2)))
var_qc_plot2


# variant mean depth
var_depth <- read_delim("variant_random_subset.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
# var_depth <- read_delim("lantana_filtered.ldepth.mean", delim = "\t",
#                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_d_stats <- summary(var_depth$mean_depth)

var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins=100)+
  ylab("count")+
  theme_few()+xlim(1,100)+
  geom_vline(xintercept = 10, color="red", linetype="dotted")+ # 10x is a good lower cutoff
  geom_vline(xintercept = 50, color="red", linetype="dotted") # upper coverage is set to remove repetitive sections or mismapping, usually uses 2x mean but this seems conservative


var_depth_plot2 <- var_depth_plot +custom_theme+
  annotate("text", size=2, color="red", x=Inf, y=Inf, hjust=1, vjust=1, label = paste("min:", round(var_d_stats[1],2),
                                                   "\nQ1:", round(var_d_stats[2],2),
                                                   "\nmedian:", round(var_d_stats[3],2),
                                                   "\nmean:", round(var_d_stats[4],2),
                                                   "\nQ3:", round(var_d_stats[5],2),
                                                   "\nmax:", round(var_d_stats[6],2)))
var_depth_plot2


# individual depth 
ind_depth <- read_delim("variant_random_subset.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# ind_depth <- read_delim("lantana_filtered.idepth", delim = "\t",
#                         col_names = c("ind", "nsites", "depth"), skip = 1)

ind_d_stats <- summary(ind_depth$depth)

ind_depth_plot <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins=100) + theme_few()+
  ylab("count")


ind_depth_plot2 <- ind_depth_plot+custom_theme+
  annotate("text", size=2,color="red", x=Inf, y=Inf, hjust=1, vjust=1,label = paste("min:", round(ind_d_stats[1],2),
                                                                         "\nQ1:", round(ind_d_stats[2],2),
                                                                         "\nmedian:", round(ind_d_stats[3],2),
                                                                         "\nmean:", round(ind_d_stats[4],2),
                                                                         "\nQ3:", round(ind_d_stats[5],2),
                                                                         "\nmax:", round(ind_d_stats[6],2)))

ind_depth_plot2
# combine plots 
comb_qc_plots <- ggarrange(var_qc_plot2, var_depth_plot2, ind_depth_plot2, nrow=3)
comb_qc_plots



# ggsave("LantCama_filtered_vcf_qc.png", comb_qc_plots, units="cm", width=12, height=12, dpi=300)

ggsave("LantCama_vcf_qc.png", comb_qc_plots, units="cm", width=12, height=12, dpi=300)
#

raw_eacp <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/eacpsamples_rawdatalookup.txt", sep="\t")
length(which((raw_eacp$targetid)!="#N/A"))

clipr::write_clip(raw_eacp$targetid)
# eacp_list <- paste(raw_eacp$targetid[which(raw_eacp$sp == "eacp")], collapse = ",")



##### Bam readcount 
bamreads <- read.csv("/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/read_counts.txt", sep="\t")
bamreads$targetid <- substr(bamreads$BAM_File, 1, 7)
bamreads$percent_unmapped <- bamreads$Unmapped_Reads/(bamreads$Unmapped_Reads + bamreads$Mapped_Reads)

bamreads2 <- merge(bamreads, meta, by="targetid")

plot_unmapped_reads <- ggplot(data=bamreads2, aes(x=cluster, y=percent_unmapped, fill=cluster))+
  geom_boxplot(outlier.size = 0.5)+theme_few()+
  scale_fill_manual(values=cluster_colours)+
  ylim(0,1)+
  labs(x="Cluster", y="Proportion unmapped reads")+
  theme(legend.position="none")

plot_unmapped_reads

ggsave("/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf_qc/bam_reads_mapping.png", plot_unmapped_reads, units="cm", width=15, height=10, dpi=300)
# #########################
# # empty because it didnt work for not diploid
# 
# # variant missingness
# var_miss <- read_delim("variant_random_subset.lmiss", delim = "\t",
#                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
# a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
# a + theme_few()
# 
# 
# # MAf
# var_freq <- read_delim("variant_random_subset.frq", delim = "\t",
#                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
# a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
# a + theme_light()
# 
# # individual missingness 
# ind_miss  <- read_delim("variant_random_subset.imiss", delim = "\t",
#                         col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
# a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
# a + theme_light()
# 
# # Ho and FIS per individual 
# ind_het <- read_delim("variant_random_subse.het", delim = "\t",
#                       col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
# a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
# a + theme_light()
