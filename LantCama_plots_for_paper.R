
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
      out_data[[paste0(species[i])]] <- data.frame(c3_species)
    }else{
      c3_species <- c3_species[!is.na(c3_species)] # remove empty columns
      out_data[[paste0(species[i])]] <- as.vector(c3_species)
      
    }
    out_data[[paste0(species[i])]] <- c3_species
}
  return(out_data)
}


test <- read_histogram_function2(m2, sp, counts2, 10) #needs meta, analysis column, counts data, and minimum number of reads per cell

####

whole_sp_plots <- function(data, species, max){
  plots <- list()
  for(i in seq_along(species)){
    x <- na.omit(unlist(list(data[[species[i]]])))
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 100, fill = "red", color="black") +
      scale_x_continuous(limits = c(0, 1),expand=c(0,0),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c(0, 1/4, 0.33, 1/2, 0.66, 3/4, 1)) +
      labs(x = "Allele reads/total locus reads", y = "", title=paste0(species[i]))+
      theme_few()+
      scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      theme(panel.grid.major.x = element_line(colour = "gray"))
    
    plots[[i]] <- p
  }
  return(plots)
}

z <- whole_sp_plots(test, c("eacp", "eawt", "per1"), c(8500,6000,2500))
ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3)


###
specific_sample_plots <- function(data, samples, max){
  plots <- list()
  for(i in seq_along(samples)){
    x <- data[samples[i],]
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 50, fill = "gray", color="black") +
      scale_x_continuous(limits = c(0, 1),expand=c(0,0),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c(0, 1/4, 0.33, 1/2, 0.66, 3/4, 1)) +
      labs(x = "Allele reads/total locus reads", y = "", title=paste0(samples[i]))+
      theme_few()+
      scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      theme(panel.grid.major.x = element_line(colour = "gray"))
    
    plots[[paste0(samples[i])]] <- p
  }
  return(plots)
}

eacp_samples <- specific_sample_plots(test$eacp, c("NSW1089413","NSW1095157","NSW1095152"), c(160,160,160))
ggarrange(eacp_samples[[1]],eacp_samples[[2]],eacp_samples[[3]], align="hv", ncol=3)


eawt_samples <- specific_sample_plots(test$eawt, c("NSW1084671","NSW1084666","NSW1095126"), c(200,200,200))
ggarrange(eawt_samples[[1]],eawt_samples[[2]],eawt_samples[[3]], align="hv", ncol=3)

eacp_samples <- specific_sample_plots(test$eacp, c("NSW1089413","NSW1095157","NSW1095152"), c(160,160,160))

