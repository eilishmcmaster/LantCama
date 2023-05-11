
read_histogram_function2 <- function(meta, sp, counts, filter_reads){
  species <- unique(meta$sp[!is.na(meta$sp)])
  
  # filter the reads
  combined_reads <- counts$c1 + counts$c2
  counts$c1[combined_reads < filter_reads] <- NA
  counts$c2[combined_reads < filter_reads] <- NA
  
  # get the proportions for all (rows are samples)
  c3_min <- pmin(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_max <- pmax(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3 <- cbind(c3_min, c3_max) 
 
  out_data  <- list() # place to put the data 
  
  for (i in seq_along(species)) {
    print(paste("Running", species[i], "now"))
    samples <- meta$sample[meta$sp == species[i]] # get the NSW ID for that species samples
    c3_species <- c3[row.names(c3) %in% samples, ] # get the readcount df with those samples
    
    if(class(c3_species) %in% c("array", "matrix")){
      c3_species <- c3_species[, colSums(!is.na(c3_species) & c3_species != "") > 0] 
      out_data[[paste0(species[i])]] <- c3_species
    }else{
      c3_species <- c3_species[!is.na(c3_species)] # remove empty columns
      out_data[[paste0(species[i])]] <- as.vector(c3_species)
      
    }
    
    out_data[[paste0(species[i])]] <- c3_species
    # par(mfrow = c(4, 5), mai = c(0.5, 0.2, 0.2, 0.2))  # Set up a 2 x 2 plotting space
    # 
    # hist(c3_species, main = species[i], xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n') # make the histogram for the species
    # axis(side = 1, at = c(0, 0.25,  0.5,  0.75, 1), labels = c(0, 0.25,  0.5,  0.75, 1))
    # 
    # if(class(c3_species) %in% c("array", "matrix")){
    #   loop.vector <- 1:nrow(c3_species)
    #   for (i in loop.vector) { # Loop over loop.vector
    #     
    #     # store data in row.i as x
    #     x <- c3_species[i,]
    #     if(sum(x, na.rm=TRUE)>0){ # skip empties
    #       # Plot histogram of x
    #       hist(x,breaks=50,
    #            main = paste(rownames(c3_species)[i]),
    #            xlab = "",#"MAF reads/ total reads",
    #            ylab="",
    #            xlim = c(0, 1),
    #            xaxt='n')
    #       axis(side = 1, at = c(0, 0.25,  0.5,  0.75, 1), labels=c(0, 0.25,  0.5,  0.75, 1))
    #     }
    #   }
    # }
}
  return(out_data)
}


test <- read_histogram_function2(m2, sp, counts2, 10) #needs meta, analysis column, counts data, and minimum number of reads per cell


library(ggplot2)

par(mfrow = c(1, 3), mai = c(0.5, 0.5, 0.2, 0.2))  # Set up a 2 x 2 plotting space

hist(test$eacp, main = "EACP", xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n') # make the histogram for the species
axis(side = 1, at = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c(0, 1/4, 0.33, 1/2, 0.66, 3/4, 1))
box()

hist(test$eawt, main = "EAWT", xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n') # make the histogram for the species
axis(side = 1, at = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c(0, 1/4, 0.33, 1/2, 0.66, 3/4, 1))
box()

hist(test$per1, main = "PER1", xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n') # make the histogram for the species
axis(side = 1, at = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c(0, 1/4, 0.33, 1/2, 0.66, 3/4, 1))
box()

if(class(c3_species) %in% c("array", "matrix")){
  loop.vector <- 1:nrow(c3_species)
  for (i in loop.vector) { # Loop over loop.vector

    # store data in row.i as x
    x <- c3_species[i,]
    if(sum(x, na.rm=TRUE)>0){ # skip empties
      # Plot histogram of x
      hist(x,breaks=50,
           main = paste(rownames(c3_species)[i]),
           xlab = "",#"MAF reads/ total reads",
           ylab="",
           xlim = c(0, 1),
           xaxt='n')
      axis(side = 1, at = c(0, 0.25,  0.5,  0.75, 1), labels=c(0, 0.25,  0.5,  0.75, 1))
  
