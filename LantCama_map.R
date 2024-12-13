# install.packages("galah")
library(galah)
library(ggrepel)
library(ggmap)
library(openxlsx)

morphid_colours <- c(pink="#EE6677", PER="forestgreen", red="red3", white="#66CCEE", orange="orange2", undetermined="#2B2B2B")
custom_theme <- theme(axis.text = element_text(size=8),
                      axis.title = element_text(size=10),
                      legend.text = element_text(size=8),
                      legend.title = element_text(size=10),
                      plot.title = element_text(size = 10),
                      legend.key.size = unit(0.5, 'lines'),
                      legend.key.height = unit(0, 'lines'))

meta <- read.xlsx("/Users/eilishmcmaster/Documents/LantCama/LantCama/outputs/LantCama_tsne_HDBSCAN_clusters.xlsx") # get meta for lat long limits

AVH_resources <- c("NSW BioNet Atlas", "PlantBank Records",
                   "NSW AVH feed", "BRI AVH data","ERBG AVH data",
                   "Centre for Australian National Biodiversity Research (CANB) AVH data",
                   "National Herbarium of Victoria (MEL) AVH data",
                   "State Herbarium of South Australia (AD) AVH data",
                   "N.C.W. Beadle Herbarium (NE) AVH data",
                   "John T. Waterhouse Herbarium (UNSW) AVH data",
                   "La Trobe University Herbarium (LTB) AVH data",
                   "Janet Cosh Herbarium (WOLL) AVH data",
                   "James Cook University Herbarium (JCT) AVH data",
                   "Tasmanian Herbarium (HO) AVH data",
                   "The University of Melbourne Herbarium (MELU) AVH data",
                   "Western Australian Herbarium (PERTH) AVH data",
                   "Northern Territory Herbarium (DNA) AVH data", "AD AVH")

galah_config(email = "eilish.summer@gmail.com") # your email here ##


fields <- c("year","country","stateProvince","decade","species",
            "stateConservation","countryConservation",
            "previousIdentifications", "associatedSequences", "basisOfRecord")

result <- galah_call() |>
  galah_identify("Lantana camara") |> # list of species to get records 
  galah_filter() |>
  galah_apply_profile(AVH) |> # list of species to get records
  galah_select(fields, group = "basic") |>
  atlas_occurrences()

result2 <- result[which(result$year >=2000 & 
                          result$country =="Australia"), ] #& result$country=="Australia"

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

whole_aus_sample_map <- ggplot()+
    geom_point(result2, 
               mapping=aes(x=decimalLongitude, y=decimalLatitude), 
               alpha=0.2, size=1, color='grey80', stroke=0)+ 
  geom_point(data=meta, 
             mapping=aes(x=long, y=lat), 
             alpha=1, size=1, color='blue', shape=16)+
    geom_sf(data=difference_poly,fill='#cbe6ef',color='transparent'
            ) +
    geom_sf(data=ozmaps::abs_ste, fill="transparent", color='black')+
    annotation_custom(
      grob = rectGrob(gp = gpar(col = "red", fill = NA, lwd = 1 )), 
      xmin = xlims3[1], xmax = xlims3[2], ymin = ylims3[1], ymax = ylims3[2]
    )+
    coord_sf(xlim=xlims2, ylim=ylims2)+
    scale_x_continuous(labels = scales::number_format(accuracy = 5),
                       breaks= seq(from = floor(min(xlims2) / 5) * 5,
                                   to = ceiling(max(xlims2) / 5) * 5,
                                   by = 5)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 5))+
  theme(panel.background = element_rect(fill="white"))+
  labs(x=element_blank(), y=element_blank())+
    theme_few()
    

whole_aus_sample_map

#### morphotype map ####

meta$morphid2 <- factor(meta$morphid2, levels=c("pink","PER", "red", "orange",  "white","undetermined"))


# Custom labeller function to capitalize the first letter of each facet label, except "PER"
capitalize_first <- function(x) {
  x <- ifelse(x == "PER", x, str_to_title(x))  # Capitalize other labels but leave "PER" as it is
  return(x)
}

morphotype_sample_map <- ggplot(difference_poly) + 
  geom_point(result2, 
             mapping = aes(x = decimalLongitude, y = decimalLatitude), 
             alpha = 0.2, size = 1, color = 'grey80', stroke = 0) + 
  geom_sf(fill = '#cbe6ef', color = 'transparent') + ##cbe6ef
  geom_sf(data = ozmaps::abs_ste, fill = "transparent", color = 'black') +
  labs(y = element_blank(), x = element_blank()) +
  theme_few() +
  theme( panel.border = element_rect()) +#axis.text.x = element_text(angle = 90),
  facet_wrap(~morphid2, labeller = labeller(morphid2 = capitalize_first)) +  # Apply the custom labeller
  scale_fill_manual(values = morphid_colours)+
  geom_point(meta, 
             mapping = aes(x = long, y = lat, fill = morphid2), 
             alpha = 1, size = 2, show.legend = F, stroke = 0.2, shape = 21, color = "white") +
  scale_x_continuous(limits = xlims3, labels = scales::number_format(accuracy = 5),
                     breaks= seq(from = floor(min(xlims3) / 5) * 5, 
                                 to = ceiling(max(xlims3) / 5) * 5, 
                                 by = 5)) +
  scale_y_continuous(limits=ylims3, labels = scales::number_format(accuracy = 5))

#### images plot #################################
# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(grid)
library(magick)

# Define paths for each image
image_paths <- list(
  Pink = "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/photos/Pink.jpg",
  PER = "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/photos/PER.jpg",
  Red = "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/photos/Red.jpg",
  Orange = "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/photos/Orange.jpg",
  White = "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/photos/White.jpg"
)

# Function to create the image plot with the specified title
create_image_plot <- function(image_path, title) {
  # Read and convert the image
  img <- image_read(image_path)
  img_grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  
  # Create ggplot with image and title
  ggplot() + 
    annotation_custom(img_grob) + 
    labs(title = title) +   # Add title
    theme_void() + 
    theme(
      plot.margin = margin(0, 1, 1, 1),           # Margins: top, right, bottom, left
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain", vjust = 1) # Preserve title capitalization
    )
}

# Create plots dynamically using file names as titles (preserving capitalization)
plots <- lapply(names(image_paths), function(name) {
  create_image_plot(image_paths[[name]], title = name)  # Use the exact name for the title
})


image_grid2 <- ggarrange(
  plotlist = plots, 
  ncol = 1, nrow = 5  # One column, five rows
)

aus_and_facet <- ggarrange(whole_aus_sample_map+theme(plot.margin = margin(5,0,-5,0))+custom_theme,
                           morphotype_sample_map+theme(plot.margin = margin(0,5,-5,0))+custom_theme,
                          nrow=2, heights=c(1,1.8))

combined_image_map_plot2 <- ggarrange(
  image_grid2 +theme(plot.margin = margin(-2,-2,2,5)),  # Margins: top, right, bottom, left
  aus_and_facet, # The main ggplot
  ncol = 2,             # Two columns
  widths = c(0.9, 3)
)

ggsave('LantCama/outputs/combined_image_map_plot.png', combined_image_map_plot2, dpi = 300, width = 14, height = 18, units = "cm", bg = "white")

#### tsne cluster map ####
meta$cluster[which(is.na(meta$cluster))] <- "Unclustered"
meta$cluster <- factor(meta$cluster, levels=c("A","B",'C','D', 'E', 'F', 'G', 'H', 'Unclustered'))

tsne_cols['Unclustered'] <- 'grey50'

cluster_sample_map <- ggplot(difference_poly) + 
  geom_sf(fill='#cbe6ef',color='transparent') +
  geom_sf(data=ozmaps::abs_ste,fill="transparent",color='black') +
  labs(y=element_blank(), x=element_blank())+
  theme_few()+
  theme(panel.border = element_rect())+#axis.text.x = element_text(angle=90), 
  geom_point(meta, 
             mapping=aes(x=long, y=lat, fill=cluster), 
             alpha=1, size=2, show.legend = F, stroke=0.1, shape=21, color="white")+
  facet_wrap(~cluster, nrow=2)+
  scale_fill_manual(values=tsne_cols) +
  scale_x_continuous(limits = xlims3, labels = scales::number_format(accuracy = 5),
                     breaks= seq(from = floor(min(xlims3) / 5) * 5, 
                         to = ceiling(max(xlims3) / 5) * 5, 
                         by = 5)) +
  scale_y_continuous(limits=ylims3, labels = scales::number_format(accuracy = 5))+
  custom_theme

ggsave('LantCama/outputs/cluster_sample_map.png', cluster_sample_map,dpi=600, width = 14, height = 12, units = "cm")


