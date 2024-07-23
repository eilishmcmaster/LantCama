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
dataset   <- "DLan22-7500"

basedir <- ""

outgroup_samples <- c("NSW1170317", "NSW1170262", "NSW1170252", "NSW1170334", "NSW1170258", "NSW1170316", "NSW1170253", "NSW1170330", "NSW1170314", "NSW1170332", "NSW1170315", "NSW1170331", "NSW1169348", "NSW1169353", "NSW1170250")


d1        <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)

table(outgroup_samples %in% d1$sample_names)
# qc1       <- report_dart_qc_stats(d1, basedir, species, dataset, threshold_missing_loci = 0.8)

# d2        <- exclude.samples(d1,by="file", excluded_sample_file = "LantCama/meta/Lcam4_sampfilt_missing80.txt")

meta      <- read.meta.data.full.analyses.df(d1, basedir, species, "DLan23-8067")
d2        <- dart.meta.data.merge(d1, meta) 
d3 <- remove.by.list(d2,meta$sample_names[which(!is.na(meta$analyses[,'EA_AM']))])
d4        <- remove.poor.quality.snps(d3, min_repro=0.99, max_missing=0.8)%>% remove.fixed.snps()
d5        <- sample.one.snp.per.locus.random(d4, seed=12345) 
d6 <- remove.by.maf(d5, 0.01)
dms <- remove.by.missingness(d6, 0.8)

samples_to_keep_80 <- dms$sample_names
length(samples_to_keep_80)
loci_to_keep_80 <- colnames(dms$gt)
length(loci_to_keep_80)
# 
# write.table(samples_to_keep_80,'LantCama/meta/samples_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)
# write.table(loci_to_keep_80,'LantCama/meta/loci_to_keep_80%.csv', row.names = FALSE, col.names = FALSE)

m2 <- d3$meta$analyses %>% as.data.frame

table((outgroup_samples %in% dms$sample_names))
outgroup_samples[which(outgroup_samples %in% dms$sample_names)]

#### write input for phylo trees ####

# for iqtree and mrbayes
dart2svdquartets(dms, RandRbase, species, dataset, add_pop=TRUE, pop=dms$sample_names)


install.packages('beastier')
library(beastier)
save_nexus_as_fasta('/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen/raw_SNPFilt_1SNPperClone/svdq/LantCama_DLan22-7500.nex',
                    '/Users/eilishmcmaster/Documents/LantCama/LantCama/popgen/raw_SNPFilt_1SNPperClone/svdq/LantCama_DLan22-7500.fasta')
