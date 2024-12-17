library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(tibble)


# nd_lea <- dart2lea(dms, RandRbase, species, dataset)
# kvalrange <- 1:20
# snmf1 <- snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 3, project = "new", CPU=8)
# 
# save(snmf1, file='LantCama/popgen/LantCama_EA_only_snmf.RData')

load(file='LantCama/popgen/LantCama_EA_only_snmf.RData')


# Function to generate a plot for a given K value
generate_plot_for_K <- function(K_value, snmf1, dms, hdb_df2) {
  best <- which.min(cross.entropy(snmf1, K = K_value))
  
  # Extract Q matrix and prepare the data frame
  qmatrix_df <- as_tibble(Q(snmf1, K = K_value, run = best)) %>%
    mutate(sample = dms$sample_names) %>%
    pivot_longer(-sample, names_to = "lea_cluster", values_to = "proportion")
  
  qmatrix_df2 <- merge(qmatrix_df, hdb_df2, by = 'sample')
  qmatrix_df2$cluster[is.na(qmatrix_df2$cluster)] <- "Unclustered"
  
  qmatrix_df2 <- qmatrix_df2 %>%
    mutate(lea_cluster = gsub("V", "", lea_cluster)) %>%
    arrange(lat) %>%
    mutate(sample = factor(sample, levels = unique(sample)))
  
  # Generate the plot
  plot <- ggplot(qmatrix_df2, aes(x = sample, y = proportion, fill = factor(lea_cluster))) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~cluster, scales = "free_x", space = "free_x") +
    theme_few() +
    scale_fill_brewer(palette = "Set3") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 11),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    ) +
    scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
    labs(x = "", y = "Ancestry Proportion", fill = "")+
    ggtitle(paste("K =", K_value))  # Add the K value as a title
  
  return(plot)
}

# List of K values to generate plots for
K_values <- c(3, 5, 7, 9, 11)  # Add more values as needed

# Generate a list of plots for each K value
plots <- lapply(K_values, function(K) generate_plot_for_K(K, snmf1, dms, hdb_df2))

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist = plots, ncol = 1, nrow = length(K_values))

# Display the combined plot
print(combined_plot)

ggsave('LantCama/outputs/SupFig_combined_lea.png', combined_plot, width = 20, height = 25,dpi=600, units = "cm")
