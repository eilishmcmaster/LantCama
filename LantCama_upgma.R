library(ggtree)
library(ggplot2)
library(ape)
library(tidytree)


upgma <- read.tree('/Users/eilishmcmaster/Documents/LantCama/upgma/treeUPGMA.tree')

upgma$node.label <- round(as.numeric(upgma$node.label),0)

ggtree_obj <- ggtree(upgma)
m2 <- read.xlsx('/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/LantCama_DLan23-8067_meta.xlsx')

x1 <- m2[m2$sample %in% ggtree_obj$data$label,]
rownames(x1) <- x1$sample

# svdq_pop_colours <-   named_list_maker(x1$svdq_pop_label, 'Spectral',11)
svdq_pop_colours <- colorRampPalette(c("grey20", "grey90"))(7)
names(svdq_pop_colours) <- 1:7
svdq_pop_colours <- c(c(eacp="#AA3377", per1="#228833",  eawt="#66CCEE"),svdq_pop_colours)
nation_colours <- named_list_maker(x1$national2, 'Paired',11)
morphid_colours <- c(pink="#AA3377", PER="#228833", red="#EE6677", white="#66CCEE", orange="#CCBB44", undetermined="#2B2B2B")


#### full tree ####

ggtree_obj <- ggtree(upgma, size=0.3) %<+% x1 

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop_label','morphid2','national2')]),
                offset=0.0006, width=.05, font.size=0,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("SVDq cluster","Country"), hjust=0) +
  scale_fill_manual(values=c(svdq_pop_colours,nation_colours, morphid_colours), na.value = "white") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size=0.8) +
  geom_rootedge(0.0005, size=0.3)+
  # geom_label2(data=upgma, aes(label=node),color='red', nudge_x = 0.0001, label.size=0, fill="transparent", size=1)+
  # geom_label2(data=upgma, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
  #             color='red', nudge_x = 0.0004, label.size=0, fill="transparent", size=2) +
  theme(legend.position = "none",plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)


custom_legend_theme <-   theme(legend.key.size = unit(0.5, 'lines'),
                        legend.title = element_text(size = 8),
                        legend.text = element_text(size = 6))
# Create separate ggplots for each fill scheme
plot_svdq_pop <- ggplot(m2, aes(x=long,y=lat, fill=svdq_pop_label)) +
  geom_tile() +
  scale_fill_manual(name = "SVDq cluster", values = svdq_pop_colours, na.value = "white") +
  custom_legend_theme

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = nation_colours, na.value = "grey90")+
  custom_legend_theme


plot_morphid2 <- ggplot(m2, aes(x=long,y=lat, fill=morphid2)) +
  geom_tile() +
  scale_fill_manual(name = "Morphotype", values = morphid_colours, na.value = "grey90")+
  custom_legend_theme

plot_prob <- ggplot(data.frame(label=upgma$node.label, x=1:length(upgma$node.label), y=1:length(upgma$node.label)),
                    aes(x=x, y=y,colour=label)) +
  geom_point() +scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)+labs(color="Probability")+custom_legend_theme


# Extract legends from the ggplots
legend_svdq_pop <- cowplot::get_legend(plot_svdq_pop)
legend_national2 <- cowplot::get_legend(plot_national2)
legend_prob <- cowplot::get_legend(plot_prob)
legend_morphid2 <- cowplot::get_legend(plot_morphid2)



legends <- cowplot::plot_grid(legend_prob, legend_morphid2, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 5/1) # Adjust aspect 
legends2 <- cowplot::plot_grid(legend_morphid2, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 3/1) # Adjust aspect 


combined_plot <- cowplot::plot_grid(hmt, legends,nrow = 1, rel_widths = c(1, 0.15))


ggsave("LantCama/outputs/LantCama_upgma_maf2.pdf",
       combined_plot, width = 25, height = 40, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_upgma_maf2.png",
       combined_plot, width = 25, height = 40, units = "cm", dpi=300)

#### collapsed tree ####

# Create nodes.identity as a named vector with names and values separately
# nodes.identity <- c("nqor", "bnor", "ghwy", "bjre", "per1_B", "per1_C", "per1_A", "red", "eacp", "eawt")
nodes.identity <- c("7", "6", "5", "4", "3", "2", "per1", "1", "eacp", "eawt")

names(nodes.identity) <- c(1030, 1050, 873, 896, 904, 938, 953, 821, 554, 673)


# # Scale clade function
# scaleMyClade <- function(.p, .node) {
#   offs <- offspring(.p$data, .node)
#   scaleClade(.p, .node, 1 / (nrow(offs) / 5 - 1))
# }

# Scale clade function proportional to the number of samples
scaleMyClade <- function(.p, .node) {
  offs <- offspring(.p$data, .node)
  num_samples <- nrow(offs)
  scaling_factor <- 1 / (log10(num_samples) + 1) # Adjust the scaling factor as needed
  scaleClade(.p, .node, scaling_factor)
}


# Collapse clade function with labels and colors
collapseMyClade <- function(.p, .node) {
  label <- nodes.identity[as.character(.node)]
  fill_color <- svdq_pop_colours[label]
  .p <- collapse(.p, .node, "max",  fill = fill_color, color="black", size=0.3)
  .p <- .p + geom_cladelabel(node = .node,align = TRUE, vjust=-0.02,offset.text=-0.0005, label = label,
                             fontsize = 2)
  return(.p)
}

nodes.to.collapse <- c(1030, 1050, 873, 896, 904, 938, 953, 821, 554, 673)

# Apply scaling and collapsing to the nodes
ggtree_obj2 <- Reduce(scaleMyClade, nodes.to.collapse, ggtree_obj)
ggtree_obj2 <- Reduce(collapseMyClade, nodes.to.collapse, ggtree_obj2)

# Add gheatmap and additional layers as required
hmt2 <- gheatmap(ggtree_obj2, as.matrix(x1[,c( 'morphid2','national2')]),#'svdq_pop_label',
                 offset = 0.001, width = .1, font.size = 0,
                 colnames_angle = 90, colnames_position = "top",
                 custom_column_labels = c("SVDq cluster", "Country"), hjust = 0) +
  scale_fill_manual(values = c(svdq_pop_colours,nation_colours, morphid_colours), na.value = "white") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size = 0.8) +
  geom_rootedge(0.0005, size = 0.3) +
  theme(legend.position = "none", plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
              color='red', nudge_x = -0.0003,nudge_y=1.1, label.size=0, fill="transparent", size=1.5)
  # geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)

# Plot the final tree
hmt2

combined_plot <- cowplot::plot_grid(hmt2, legends2,nrow = 1, rel_widths = c(1, 0.2))


ggsave("LantCama/outputs/LantCama_upgma_maf2_collapsed.pdf",
       combined_plot, width = 15, height = 20, units = "cm", dpi=600)
 

ggsave("LantCama/outputs/LantCama_upgma_maf2_collapsed.png",
       combined_plot, width = 15, height = 20, units = "cm", dpi=600)

#### networ ####

Nnet <- phangorn::read.nexus.networx('LantCama/outputs/all_LantCama_nexus_file_for_R_80missing_maf2.nex')


x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)
x$svdq_pop_label[!is.na(x$sample)&is.na(x$svdq_pop_label)] <- "ungrouped"

net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)


levels2 <- unique(dms$meta$analyses[,"national2"]) %>% sort()
levels2_shapes <- setNames(1:nlevels(factor(levels2)), levels2)

hull <- x %>% group_by(svdq_pop_label) %>% 
  slice(chull(x, y))

hull <- hull[!hull$svdq_pop_label=='ungrouped',]

x_filtered <- x[!is.na(x$national2),]

splitstree_plot_svdq <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.7, expand = 0.01, radius = 0.01,
             aes(fill = svdq_pop_label, color = "transparent")) +
  geom_point(data = x_filtered, aes(shape = national2), color="white", size=2) +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_filtered, aes(color = morphid2, shape = national2)) +
  scale_fill_manual(values = svdq_pop_colours, na.translate = FALSE) +
  scale_colour_manual(values = morphid_colours, na.translate = FALSE) +
  theme_void() +
  expand_limits(x = c(min(x_filtered$x) - 0.01 * net_x_axis, max(x_filtered$x) + 0.01 * net_x_axis),
                y = c(min(x_filtered$y) - 0.01 * net_y_axis, max(x_filtered$y) + 0.01 * net_y_axis)) +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  coord_fixed() +
  labs(color = "Morphotype", shape = "Origin", fill = "SVDq clusters") +
  guides(colour = guide_legend(title.position = "top", nrow = 2, override.aes = list(fill = NA, linetype = 0)), 
         fill = guide_legend(title.position = "top", nrow = 2, override.aes = list(color = NA)),
         shape = guide_legend(title.position = "top", nrow = 2))
# Save the plot
ggsave("LantCama/outputs/LantCama_splitstree_svdqpops_80miss_maf2.pdf",
       splitstree_plot_svdq, width = 20, height = 20, units = "cm", dpi = 600)

ggsave("LantCama/outputs/LantCama_splitstree_svdqpops_80miss_maf2.png",
       splitstree_plot_svdq, width = 20, height = 20, units = "cm", dpi = 600)

