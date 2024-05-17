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

svdq_pop_colours <-   named_list_maker(x1$svdq_pop, 'Spectral',11)
nation_colours <- named_list_maker(x1$national2, 'Paired',11)


#### full tree ####

ggtree_obj <- ggtree(upgma, size=0.3, layout="roundrect") %<+% x1 

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop','national2')]),
                offset=0.0006, width=.05, font.size=0,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("SVDq cluster","Country"), hjust=0) +
  scale_fill_manual(values=c(svdq_pop_colours,nation_colours), na.value = "grey90") +
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
plot_svdq_pop <- ggplot(m2, aes(x=long,y=lat, fill=svdq_pop)) +
  geom_tile() +
  scale_fill_manual(name = "SVDq cluster", values = svdq_pop_colours, na.value = "grey90") +
  custom_legend_theme

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = nation_colours, na.value = "grey90")+
  custom_legend_theme

plot_prob <- ggplot(data.frame(label=upgma$node.label, x=1:length(upgma$node.label), y=1:length(upgma$node.label)),
                    aes(x=x, y=y,colour=label)) +
  geom_point() +scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)+labs(color="Probability")+custom_legend_theme


# Extract legends from the ggplots
legend_svdq_pop <- cowplot::get_legend(plot_svdq_pop)
legend_national2 <- cowplot::get_legend(plot_national2)
legend_prob <- cowplot::get_legend(plot_prob)

legends <- cowplot::plot_grid(legend_prob,legend_svdq_pop, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 6/1) # Adjust aspect 
legends2 <- cowplot::plot_grid(legend_svdq_pop, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 2/1) # Adjust aspect 


combined_plot <- cowplot::plot_grid(hmt, legends,nrow = 1, rel_widths = c(1, 0.15))


ggsave("LantCama/outputs/LantCama_upgma_maf2.pdf",
       combined_plot, width = 25, height = 40, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_upgma_maf2.png",
       combined_plot, width = 25, height = 40, units = "cm", dpi=300)

#### collapsed tree ####

# Create nodes.identity as a named vector with names and values separately
nodes.identity <- c("nqor", "bnor", "ghwy", "bjre", "per1_B", "per1_C", "per1_A", "red", "eacp", "eawt")
names(nodes.identity) <- c(1030, 1050, 873, 896, 904, 938, 953, 821, 554, 673)


# Scale clade function
scaleMyClade <- function(.p, .node) {
  offs <- offspring(.p$data, .node)
  scaleClade(.p, .node, 1 / (nrow(offs) / 3 - 1))
}

# Collapse clade function with labels and colors
collapseMyClade <- function(.p, .node) {
  label <- nodes.identity[as.character(.node)]
  fill_color <- svdq_pop_colours[label]
  .p <- collapse(.p, .node, "max",  fill = fill_color, color="black", size=0.3)
  # .p <- .p + geom_cladelabel(node = .node,align = TRUE, vjust=-0.02, label = label, 
  #                            fontsize = 2)
  return(.p)
}

nodes.to.collapse <- c(1030, 1050, 873, 896, 904, 938, 953, 821, 554, 673)

# Apply scaling and collapsing to the nodes
ggtree_obj2 <- Reduce(scaleMyClade, nodes.to.collapse, ggtree_obj)
ggtree_obj2 <- Reduce(collapseMyClade, nodes.to.collapse, ggtree_obj2)

# Add gheatmap and additional layers as required
hmt2 <- gheatmap(ggtree_obj2, as.matrix(x1[,c('svdq_pop','national2')]),
                 offset = 0.001, width = .05, font.size = 0,
                 colnames_angle = 90, colnames_position = "top",
                 custom_column_labels = c("SVDq cluster", "Country"), hjust = 0) +
  scale_fill_manual(values = c(cc, cc2), na.value = "grey90") +
  theme_tree2() +
  geom_tiplab(aes(label = label), size = 0.8) +
  geom_rootedge(0.0005, size = 0.3) +
  theme(legend.position = "none", plot.margin = margin(0, -0.8, 0, 0, "cm"), axis.text.x = element_text(size=6))+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80),
              color='red', nudge_x = -0.0003,nudge_y=1.1, label.size=0, fill="transparent", size=1.5)
  # geom_nodepoint(aes(color=as.numeric(label)), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)

# Plot the final tree
hmt2

combined_plot <- cowplot::plot_grid(hmt2, legends2,nrow = 1, rel_widths = c(1, 0.2))


ggsave("LantCama/outputs/LantCama_upgma_maf2_collapsed.pdf",
       combined_plot, width = 15, height = 20, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_upgma_maf2_collapsed.png",
       combined_plot, width = 15, height = 20, units = "cm", dpi=600)
