library(dplyr)
library(ggpubr)
library(popkin)
library(SNPRelate)
library(vcfR)

### Import data ###

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/kin_compare_functions.R?raw=TRUE")

vcfraw <- read.vcfR('/Users/eilishmcmaster/Documents/LantCama/LantCama/vcf/raw_incomp.mac3bcf.vcf')
# genind <- vcfR2genind(vcf)
# gt <- genind$tab

ad <- extract.gt(vcfraw, element = 'AD')

#remove homozygous loci
gt <- extract.gt(vcfraw, element = 'GT')
hets <- is_het(gt)

is.na( ad[ !hets ] ) <- TRUE


df <- ad

