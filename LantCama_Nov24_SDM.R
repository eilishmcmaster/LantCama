library(ggplot2)
library(sf)
library(ggthemes)
library(terra)
library(openxlsx)
library(grid)
library(ggpubr)

custom_theme <- theme(axis.text = element_text(size=8),
                      axis.title = element_text(size=10),
                      legend.text = element_text(size=8),
                      legend.title = element_text(size=10),
                      plot.title = element_text(size = 10),
                      legend.key.size = unit(0.5, 'lines'),
                      legend.key.height = unit(0.5, 'lines'))
#### ausmap ####

### invert australia ####
# Define the bounding box as a polygon
xlims <- c(96.81703, 167.9969)
ylims <- c(-8, -45)

xlims2 <- c(112, 155)
ylims2 <- c(-10, -43)

### bound for zoom 

xlims3 <- c(140, 155)
ylims3 <- c(-36.4,-14.5)

#### inverse australia polygon ####
# Create the bounding box polygon
bbox_poly <- st_as_sf(st_sfc(st_polygon(list(matrix(c(
  xlims[1], ylims[1],
  xlims[2], ylims[1],
  xlims[2], ylims[2],
  xlims[1], ylims[2],
  xlims[1], ylims[1]
), ncol = 2, byrow = TRUE))), crs = 4326)) # Use CRS consistent with GDA94

# Load ozmaps::abs_ste
ozmaps::abs_ste -> australia_map

# Ensure consistent CRS
australia_map <- st_transform(australia_map, st_crs(bbox_poly))

# Perform the difference operation
difference_poly <- st_difference(bbox_poly, st_union(australia_map))

#### import data ####

meta <- read.xlsx('LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx')


# ggplot()+
#   geom_point(data=thomas2006, mapping=aes(x=long, y=lat, color=morphid))+
#   geom_sf(data=ozmaps::abs_ste, fill="transparent", color='black')
  
thomas2006 <- read.xlsx('LantCama/meta/Thomas2006_Lantana_flower_data.xlsx')
thomas2006_2 <- thomas2006[thomas2006$morphid %in% c('pink'),]



# Define file paths for each group
file_paths <- list(
  A = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_A/maxent_predict.grd',
  B = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_B/maxent_predict.grd',
  C = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_C/maxent_predict.grd',
  D = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_D/maxent_predict.grd',
  E = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_E/maxent_predict.grd',
  `F` = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_F/maxent_predict.grd',
  G = '/Users/eilishmcmaster/Documents/LantCama/lantana_maxent_results_v3/lantana_v3_G/maxent_predict.grd'
)

# Load and process rasters for each group
data_list <- lapply(file_paths, function(path) {
  rast_data <- rast(path)                             # Load the raster
  df <- as.data.frame(rast_data, xy = TRUE)          # Convert to data frame
  df <- df[which(df$layer > 0.01), ]                 # Filter based on layer > 0.01
  return(df)
})

# Assign each processed data frame to a variable
names(data_list) <- names(file_paths)
max_df_a <- data_list$A
max_df_b <- data_list$B
max_df_c <- data_list$C
max_df_d <- data_list$D
max_df_e <- data_list$E
max_df_f <- data_list$`F`
max_df_g <- data_list$G

# Check one of the resulting data frames (e.g., max_df_a)
head(max_df_a)


# Function to create a plot for each group
create_plot <- function(r_df, title) {
  ggplot() +
    geom_raster(data = r_df, mapping = aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "blue", high = "darkblue",
                         midpoint = 0.5, limits = c(0, 1), na.value = 'white') +
    geom_sf(data = difference_poly, fill = '#cbe6ef', color = 'transparent') +
    geom_sf(data = ozmaps::abs_ste, fill = "transparent", color = 'black') +
    coord_sf(xlim = xlims2, ylim = ylims2) +
    scale_x_continuous(labels = scales::number_format(accuracy = 5),
                       breaks = seq(from = floor(min(xlims2) / 5) * 5,
                                    to = ceiling(max(xlims2) / 5) * 5,
                                    by = 5)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 5)) +
    theme(panel.background = element_rect(fill = "white")) +
    labs(x = element_blank(), y = element_blank()) +
    theme_few() + custom_theme +
    labs(fill = "Habitat\nsuitability", title=title)
}

# thomas2006$Rating <- factor(thomas2006$Rating, levels = c("0", "1", "2", "3"))

# Create and display plots for A, B, F, and D
plot_a <- create_plot(max_df_a, "Cluster A (Common Pink)")
plot_a2 <- plot_a+
  geom_point(data = thomas2006[thomas2006$morphid=="pink",], 
             mapping = aes(x = long, y = lat), 
             shape = 20, size = 1.5) + 
  geom_point(data = thomas2006_2, 
             mapping = aes(x = long, y = lat, color = as.factor(Rating)), 
             shape = 20, size = 1) + 
  scale_color_manual(values = c(`0` = "red", `3` = "green", `2` = "yellow", `1` = "orange"),
                     drop=TRUE)+
  labs(color="Susceptibility")

plot_b <- create_plot(max_df_b, "Cluster B (Common PER)")
plot_c <- create_plot(max_df_c, "Cluster C (Oblong Red)")
plot_d <- create_plot(max_df_d, "Cluster D (Helidon White)")
plot_e <- create_plot(max_df_e, "Cluster E (True Orange)")
plot_f <- create_plot(max_df_f, "Cluster F (Townsville Red-centred Pink)")
plot_g <- create_plot(max_df_g, "Cluster G (Townsville Prickly Orange)")


# # Display the plots
# 
# combined_sdms <- ggarrange(plot_a2,
#           plot_f,
#           plot_b,
#           plot_d,
#           ncol=2,
#           # nrow=2, ncol=2, 
#           align='hv',
#           common.legend=T, legend="right", 
#           abels="auto") #font.label = list(face = "plain")
# 
# ggsave('LantCama/outputs/LantCama_combined_sdms.png', plot=(combined_sdms), 
#        dpi = 300,
#        width = 20, height = 20, units = "cm")



library(gridExtra)

# Ensure combined_sdms is created with ggarrange
combined_sdms <- ggarrange(
  plot_a2, plot_f,
  plot_b, plot_d, # Replace with your actual plots
  ncol = 2, nrow = 2,
  align='hv',
  common.legend=T, legend="right",
  labels = 'auto', font.label = list(face = "plain")
)

# Save the plot using ggsave
ggsave(
  filename = "LantCama/outputs/LantCama_combined_sdms.png",
  plot = combined_sdms,
  dpi = 300,
  width = 20,
  height = 16,
  units = "cm", bg = 'white'
)


suit_legend <- cowplot::get_legend(plot_f)

combined_sdms2 <- cowplot::plot_grid(
                  plot_a + theme(legend.position = 'none'), 
                   plot_b + theme(legend.position = 'none'),
                   plot_c + theme(legend.position = 'none'), 
                   plot_d + theme(legend.position = 'none'),
                   plot_e + theme(legend.position = 'none'), 
                   plot_f + theme(legend.position = 'none'), 
                   plot_g + theme(legend.position = 'none'),
                   suit_legend,
                   ncol=3,nrow=3,
                   align='hv',
                   # nrow=2, ncol=2, align='hv',
                   labels=c('a','b','c','d','e','f','g'),
                  label_fontface ="plain")


# Save the plot using ggsave
ggsave(
  filename = "LantCama/outputs/LantCama_combined_sdms2.png",
  plot = combined_sdms2,
  dpi = 300,
  width = 20,
  height = 20,
  units = "cm", bg = 'white'
)
