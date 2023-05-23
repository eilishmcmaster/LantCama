
read_histogram_function2 <- function(meta, counts, filter_reads, species_col) {
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
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
    samples <- meta$sample[meta[[species_col]] == species[i]] # get the NSW ID for that species samples
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


test <- read_histogram_function2(m2, counts2, 10, species_col="sp") #needs meta, analysis column, counts data, and minimum number of reads per cell

####

whole_sp_plots <- function(data, species, max){
  plots <- list()
  for(i in seq_along(species)){
    x <- na.omit(unlist(list(data[[species[i]]])))
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 100, fill = "green", color="black", size=0.2) + 
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4", "1/3", "1/2", "2/3", "3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(species[i]))+
      theme_few()+
      scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black")),
            axis.text.y = element_text(size=8),
            title = element_text(size=10),
            plot.margin = margin(3, 1, 1, 1))
    
    plots[[i]] <- p
  }
  return(plots)
}

z <- whole_sp_plots(test, c("eacp", "eawt", "per1"), c(8500,8500,8500))
ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3,
          labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")
ggsave("LantCama/outputs/combined_group.png", plot = last_plot(), width = 190, height = 80, dpi = 300, units = "mm")


###
specific_sample_plots <- function(data, samples, max){
  plots <- list()
  for(i in seq_along(samples)){
    x <- data[samples[i],]
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 30, fill = "gray", color="black", size=0.2) + #usually use 50 bins
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4 ", "1/3", "1/2", "2/3", " 3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(samples[i]))+
      theme_few()+
      scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black")),
            axis.text.y = element_text(size=8),
            title = element_text(size=8),
            plot.margin = margin(3, 1, 1, 3))

    
    plots[[paste0(samples[i])]] <- p
  }
  return(plots)
}

eacp_samples <- specific_sample_plots(test$eacp,
                                      c("NSW1089413","NSW1096776","NSW1095152"),#"NSW1095157"
                                      rep(250,3))#c(160,160,160)) # for 50 breaks

eawt_samples <- specific_sample_plots(test$eawt,
                                      c("NSW1084671","NSW1084666","NSW1095126"),
                                      rep(350,3))#c(250,250,250))

per1_samples <- specific_sample_plots(test$per1,
                                      c("NSW1158953","NSW1150367","NSW1161296"), #"NSW1152374"),
                                      rep(80,3))#c(50,50,50))

ggarrange(eacp_samples[[1]],eacp_samples[[2]],eacp_samples[[3]],
          eawt_samples[[1]],eawt_samples[[2]],eawt_samples[[3]],
          per1_samples[[1]],per1_samples[[2]],per1_samples[[3]],
          align="hv", ncol=3, nrow=3,
          labels=c("A","","","B","","","C"),
          font.label = list(size = 10, color = "black", face = "bold", family = NULL))%>%
annotate_figure(.,
  bottom = "Allele frequency",
  left="Count"
)
# ggsave("LantCama/outputs/specific_examples.png", plot = last_plot(), width = 190, height = 160, dpi = 300, units = "mm")

