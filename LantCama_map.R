# install.packages("galah")
library(galah)
library(ggrepel)
library(ggmap)
library(openxlsx)

morphid_colours <- c(pink="#EE6677", PER="forestgreen", red="red3", white="#66CCEE", orange="orange2", undetermined="#2B2B2B")


# meta <- read.xlsx("/Users/eilishmcmaster/Documents/ZierObco/ZierObco/meta/ZierObco_DZ22-7321_meta.xlsx") # get meta for lat long limits

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

divxlims <- c(min(result2$decimalLongitude, na.rm=TRUE)-0.3,
              max(result2$decimalLongitude, na.rm=TRUE)+0.3) #find the min / max longitude

divylims <- c(min(result2$decimalLatitude, na.rm=TRUE)-0.2,
              max(result2$decimalLatitude, na.rm=TRUE)+0.2) #find the min / max latitude

xlims <- c(113, 155)
ylims <- c(-43,-10)

library("ozmaps")
test_map <- ggplot(ozmaps::abs_ste) + 
  geom_sf(fill="white", color='black') +
  geom_point(result2, 
             mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             alpha=1, size=2, color='grey80')+ #steelblue
  geom_sf(alpha=0, color='black') +
  # coord_sf(xlim = divxlims, ylim = divylims) +
  xlim(xlims)+
  ylim(ylims)+
  labs(y=element_blank(), x=element_blank(), colour="ALA records")+
  theme_few()+
  theme(axis.text.x = element_text(angle=90), panel.background = element_rect(fill="grey90"))+
  geom_point(m2, 
             mapping=aes(x=long, y=lat), 
             alpha=1, size=1, color='red', shape=16)

test_map

### invert australia ####
# Define the bounding box as a polygon
xlims2 <- c(96.81703, 167.9969)
ylims2 <- c(-8, -45)

library(sf)
# Create the bounding box polygon
bbox_poly <- st_as_sf(st_sfc(st_polygon(list(matrix(c(
  xlims2[1], ylims2[1],
  xlims2[2], ylims2[1],
  xlims2[2], ylims2[2],
  xlims2[1], ylims2[2],
  xlims2[1], ylims2[1]
), ncol = 2, byrow = TRUE))), crs = 4326)) # Use CRS consistent with GDA94

# Load ozmaps::abs_ste
ozmaps::abs_ste -> australia_map

# Ensure consistent CRS
australia_map <- st_transform(australia_map, st_crs(bbox_poly))

# Perform the difference operation
difference_poly <- st_difference(bbox_poly, st_union(australia_map))

whole_aus_sample_map <- ggplot(difference_poly) + 
  # geom_sf(fill="white", color='black') +
  geom_point(result2, 
             mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             alpha=0.2, size=1, color='grey80', stroke=0)+ 
  geom_sf(fill='lightblue2',color='transparent') +
  xlim(xlims)+
  ylim(ylims)+
  labs(y=element_blank(), x=element_blank(), colour="ALA records")+
  theme_few()+
  theme(panel.background = element_rect(fill="white"))+
  geom_sf(data=ozmaps::abs_ste, fill="transparent", color='black') +
  geom_point(m2, 
             mapping=aes(x=long, y=lat), 
             alpha=1, size=1, color='red', shape=16)


#
ggsave('LantCama/outputs/test_map.pdf', whole_aus_sample_map, width = 20, height = 15, units = "cm")


xlims2 <- c(140, 155)
ylims2 <- c(-37,-10)

m2$morphid2 <- factor(m2$morphid2, levels=c("pink","PER", "red", "orange",  "white","undetermined"))


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
  xlim(xlims2) +
  ylim(ylims2) +
  labs(y = element_blank(), x = element_blank()) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90), panel.border = element_rect()) +
  geom_point(m2, 
             mapping = aes(x = long, y = lat, fill = morphid2), 
             alpha = 1, size = 2, show.legend = F, stroke = 0.1, shape = 21, color = "white") +
  facet_wrap(~morphid2, labeller = labeller(morphid2 = capitalize_first)) +  # Apply the custom labeller
  scale_fill_manual(values = morphid_colours)

# Print the plot
print(morphotype_sample_map)


# ggsave('LantCama/outputs/morphotype_sample_map.pdf', morphotype_sample_map, width = 20, height = 20, units = "cm")
ggsave('LantCama/outputs/morphotype_sample_map.png', morphotype_sample_map,dpi=600, width = 20, height = 20, units = "cm")


ggarrange(whole_aus_sample_map, morphotype_sample_map, nrow=2, heights=c(1.2,2))


morphotype_sample_map2 <- ggplot(difference_poly) + 
  geom_point(result2, 
             mapping = aes(x = decimalLongitude, y = decimalLatitude), 
             alpha = 0.2, size = 1, color = 'grey80', stroke = 0) + 
  geom_sf(fill = '#cbe6ef', color = 'transparent') + ##cbe6ef
  geom_sf(data = ozmaps::abs_ste, fill = "transparent", color = 'black') +
  xlim(xlims2) +
  ylim(ylims2) +
  labs(y = element_blank(), x = element_blank()) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90), panel.border = element_rect()) +
  geom_point(m2[which(m2$morphid2!="undetermined"),], 
             mapping = aes(x = long, y = lat, fill = morphid2), 
             alpha = 1, size = 2, show.legend = F, stroke = 0.1, shape = 21, color = "white") +
  facet_wrap(~morphid2, labeller = labeller(morphid2 = capitalize_first), nrow = 1) +  # Apply the custom labeller
  scale_fill_manual(values = morphid_colours)

# ####
# data_2003 <- read.csv('/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/Lantana distribution 2003 paper cleaned.csv')
# colnames(data_2003)[3] <- 'morphid2'
# data_2003$morphid2 <- factor(data_2003$morphid2, levels=c("pink","PER", "red", "orange",  "white","undetermined"))
# 
# ggplot(difference_poly) + 
#   geom_point(result2,
#              mapping = aes(x = decimalLongitude, y = decimalLatitude),
#              alpha = 0.2, size = 1, color = 'grey80', stroke = 0) +
#   geom_point(data_2003,
#              mapping = aes(x = long, y = lat),
#              alpha = 1, size = 1, color = 'grey20', stroke = 0) +
#   geom_sf(fill = '#cbe6ef', color = 'transparent') + ##cbe6ef
#   geom_sf(data = ozmaps::abs_ste, fill = "transparent", color = 'black') +
#   xlim(xlims2) +
#   ylim(ylims2) +
#   labs(y = element_blank(), x = element_blank()) +
#   theme_few() +
#   theme(axis.text.x = element_text(angle = 90), panel.border = element_rect()) +
#   geom_point(m2, 
#              mapping = aes(x = long, y = lat, color = morphid2), 
#              alpha = 1, size = 1, show.legend = F, shape = 4, stroke=1) +
#   facet_wrap(~morphid2, labeller = labeller(morphid2 = capitalize_first)) +  # Apply the custom labeller
#   scale_color_manual(values = morphid_colours)

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
      plot.margin = margin(0, 1, 0, 1),           # Margins: top, right, bottom, left
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain", vjust = 1) # Preserve title capitalization
    )
}

# Create plots dynamically using file names as titles (preserving capitalization)
plots <- lapply(names(image_paths), function(name) {
  create_image_plot(image_paths[[name]], title = name)  # Use the exact name for the title
})

# Arrange the plots in a single column
image_grid <- ggarrange(
  plotlist = plots, 
  nrow = 1, ncol = 5  # One column, five rows
)

# Assuming morphotype_sample_map is already defined, combine it with the image grid
combined_plot <- ggarrange(
  image_grid+theme(plot.margin = margin(0,10,0,40)),           # Image grid plot
  morphotype_sample_map2, # The main ggplot
  nrow = 2,             # Two columns
  heights = c(1, 2.5)
  )
  # labels = c("a", "b"),
  # font.label = list(face = "plain"))

# Save the combined plot to a file
ggsave('LantCama/outputs/morphotype_sample_map.png', combined_plot, dpi = 300, width = 15, height = 10, units = "cm", bg = "white")

#### tsne cluster map ####
m2$cluster <- factor(m2$cluster, levels=c("A","B",'C','D', 'E', 'F', 'G', 'H', NA))

cluster_sample_map <- ggplot(difference_poly) + 
  # geom_point(result2, 
  #            mapping=aes(x=decimalLongitude, y=decimalLatitude), 
  #            alpha=0.2, size=1, color='grey80', stroke=0)+ 
  geom_sf(fill='#cbe6ef',color='transparent') +
  geom_sf(data=ozmaps::abs_ste,fill="transparent",color='black') +
  xlim(xlims2)+
  ylim(ylims2)+
  labs(y=element_blank(), x=element_blank())+
  theme_few()+
  theme(axis.text.x = element_text(angle=90), panel.border = element_rect())+
  geom_point(m2, 
             mapping=aes(x=long, y=lat, fill=cluster), 
             alpha=1, size=2, show.legend = F, stroke=0.1, shape=21, color="white")+
  facet_wrap(~cluster, nrow=2)+
  scale_fill_manual(values=tsne_cols)

ggsave('LantCama/outputs/cluster_sample_map.png', cluster_sample_map,dpi=600, width = 20, height = 19, units = "cm")


# cluster_sample_map <- ggplot(difference_poly) + 
#   # geom_point(result2, 
#   #            mapping=aes(x=decimalLongitude, y=decimalLatitude), 
#   #            alpha=0.2, size=1, color='grey80', stroke=0)+ 
#   geom_sf(fill='#cbe6ef',color='transparent') +
#   geom_sf(data=ozmaps::abs_ste,fill="transparent",color='black') +
#   xlim(xlims2)+
#   ylim(ylims2)+
#   labs(y=element_blank(), x=element_blank())+
#   theme_few()+
#   theme(axis.text.x = element_text(angle=90), panel.border = element_rect())+
#   geom_point(m2, 
#              mapping=aes(x=long, y=lat, fill=cluster), 
#              alpha=1, size=2, show.legend = F, stroke=0.1, shape=21, color="white")+
#   # facet_wrap(~cluster, nrow=2)+
#   scale_fill_manual(values=tsne_cols2, na.value = NA)

