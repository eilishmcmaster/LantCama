
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

counts2 <- read_dart_counts_csv_faster('LantCama/dart_raw/Report-DLan22-7500/Report_DLan22_7500_7160part_6495_6370/Report_DLan22-7500_3_moreOrders_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=1, 
                                       minGenotypeCount=0)

d1 <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE, altcount = TRUE) #read DArT file

# missingness threshold is 1 to prevent alleles of small species groups being removed 
d2 <- remove.poor.quality.snps(d1, min_repro=0.96, max_missing=missingness)

d3 <- remove.fixed.snps(d2)
d3 <- sample.one.snp.per.locus.random(d3, seed=12345)

m2 <- custom.read(species, dataset) #read custom metadata csv

mm1 <- read.meta.data(d3, RandRbase, species, dataset, fields=(ncol(m2)-4))

dms       <- dart.meta.data.merge(d3, mm1)


#plot function 
a_function <- function(dms, counts){
  species <- unique(dms$meta$analyses[,"sp"])
  print(species)
  
  for (i in 1:length(species)) {
    print(paste("Running", species[i],"now"))
    tryCatch({
      dmsx <- remove.by.list(dms, m2[(m2$sp %in% paste(species[i])),] %>%.$sample) %>% remove.by.maf(., 0.05)
      test2 <- count_subsetter(dmsx, counts, 0.05)
    }, error = function(e) {
      message("Error: not enough loci")
    })
    if(length(counts$sample_names)>0){
      tr <-  t(test2$c1)
      o <- lapply(split(tr,rownames(tr)), as.list)
      
      tr2 <-  t(test2$c2)
      o2 <- lapply(split(tr2,rownames(tr2)), as.list)
      
      nn <- mergeLists_internal(o, o2)
      
      minor <- lapply(nn, sapply, function(x) min(x)/sum(x))
      a <- do.call(rbind, minor) #make matrix
      major <- lapply(nn, sapply, function(x) max(x)/sum(x))
      b <- do.call(rbind, major) #make matrix
      c <- cbind(a,b)
      
      # only for six samples at a time
      par(mfrow = c(4, 5), mai=c(0.5,0.2,0.2,0.2))  # Set up a 2 x 2 plotting space
      
      # Create the loop.vector (all the columns)
      loop.vector <- 1:nrow(c)
      z <- paste(unique(dms$meta$analyses[,"sp"])[i])
      
      hist(c, main=z, xlab="", ylab="", breaks=50, col="red")
      
      for (i in loop.vector) { # Loop over loop.vector
        
        # store data in row.i as x
        x <- c[i,]
        if(sum(x, na.rm=TRUE)>0){ # skip empties
          # Plot histogram of x
          hist(x,breaks=50,
               main = paste(rownames(c)[i]),
               xlab = "",#"MAF reads/ total reads",
               ylab="",
               xlim = c(0, 1))
        }
      }
      
      
    }
  }
  
  
}

dms_filter<- remove.by.list(dms, m2[!is.na(m2$sp),] %>%.$sample) 

counts <- count_subsetter(dms_filter, counts2, 0.00)


pdf(file="lantana_all.pdf")
a_function(dms_filter, counts)
dev.off()

