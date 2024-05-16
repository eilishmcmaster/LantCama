library(ggtree)
library(ggplot2)

mrbayes <- treeio::read.mrbayes('/Users/eilishmcmaster/Documents/LantCama/mrbayes/lantana_eilish_mrbayes_full/infile.nex.con.tre')
mrbayes@data[["prob"]] <- round(as.numeric(mrbayes@data[["prob"]]), 2)

ggtree_obj <- ggtree(mrbayes)
m2 <- read.xlsx('/Users/eilishmcmaster/Documents/LantCama/LantCama/meta/LantCama_DLan23-8067_meta.xlsx')

x1 <- m2[m2$sample %in% ggtree_obj$data$label,]
rownames(x1) <- x1$sample

cc <-   named_list_maker(x1$svdq_pop, 'Spectral',11)
cc2 <- named_list_maker(x1$national2, 'Paired',11)

ggtree_obj <- ggtree(mrbayes, size=0.3) %<+% x1 
# ggtree_obj <- ggtree_obj%>% scaleClade(545, .2) %>%
#   collapse(545,"max")

hmt <- gheatmap(ggtree_obj, as.matrix(x1[,c('svdq_pop','national2')]),
                offset=0.025, width=.05, font.size=0,
                colnames_angle=90, colnames_position="top",
                custom_column_labels=c("SVDq cluster","Country"), hjust=0) +
  scale_fill_manual(values=c(cc,cc2), na.value = "grey90") +
  theme_tree2() +
  # geom_rootedge(0.0005, size=0.01) +
  geom_tiplab(aes(label = label), size=0.8, align=TRUE, linetype='dotted', linesize=0.2) +
  # geom_label2(data=mrbayes, aes(label=prob),
  #             color='red', nudge_x = 0.00015, label.size=0, fill="transparent", size=1) +
  theme(legend.position = "none",plot.margin = margin(0, -0.8, 0, 0, "cm")) +
  # scale_y_continuous(expand=c(0.05, 0))+
  # geom_text2(data=mrbayes, aes(label=node), size=0.5, color="red")+
  geom_nodepoint(data=mrbayes, aes(color=prob), size=0.5)+scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)
# geom_hilight(node=600, fill="steelblue", alpha=0.5) +
# geom_hilight(node=624, fill="darkgreen", alpha=0.5) +
# geom_hilight(node=716, fill="gray", alpha=0.5) +
# geom_hilight(node=734, fill="pink", alpha=0.5) +
# geom_hilight(node=729, fill="beige", alpha=0.5) +
# geom_hilight(node=545, fill="yellow", alpha=0.5) +
# geom_hilight(node=562, fill="blue", alpha=0.5) 

# Create separate ggplots for each fill scheme
plot_svdq_pop <- ggplot(m2, aes(x=long,y=lat, fill=svdq_pop)) +
  geom_tile() +
  scale_fill_manual(name = "SVDq cluster", values = cc, na.value = "grey90") 

plot_national2 <- ggplot(m2, aes(x=long,y=lat, fill=national2)) +
  geom_tile() +
  scale_fill_manual(name = "Country", values = cc2, na.value = "grey90")

plot_prob <- ggplot(data.frame(prob=mrbayes@data[["prob"]], x=1:length(mrbayes@data[["prob"]]), y=1:length(mrbayes@data[["prob"]])),
                    aes(x=x, y=y,colour=prob)) +
  geom_point() +scale_color_distiller(palette = "RdYlBu", direction=+1, na.value = NA)+labs(color="Probability")

# Extract legends from the ggplots
legend_svdq_pop <- cowplot::get_legend(plot_svdq_pop)
legend_national2 <- cowplot::get_legend(plot_national2)
legend_prob <- cowplot::get_legend(plot_prob)

legends <- cowplot::plot_grid(legend_prob,legend_svdq_pop, legend_national2, ncol=1, align="hv") + theme(aspect.ratio = 6/1) # Adjust aspect 


combined_plot <- cowplot::plot_grid(hmt, legends,nrow = 1, rel_widths = c(1, 0.15))


ggsave("LantCama/outputs/LantCama_mrbayes_maf2.pdf",
       combined_plot, width = 25, height = 40, units = "cm", dpi=600)

ggsave("LantCama/outputs/LantCama_mrbayes_maf2.png",
       combined_plot, width = 25, height = 40, units = "cm", dpi=300)

