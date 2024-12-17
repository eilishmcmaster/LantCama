hybrids_only <- c(meta$sample_names[which(!is.na(meta$analyses[,'hybsites']))])

pca_hybrids <- g_pca_df2[which(g_pca_df2$Row.names %in% hybrids_only),]

ggplot(pca_hybrids,
       aes(x=PC1, y=PC2, colour=morphid2))+ xlab(pcnames[1])+ylab(pcnames[2])+
  geom_point(size=0.5)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="HBDSCAN Cluster")+
  theme(legend.key.size = unit(0, 'lines'))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)


ind_ho <- fastDiversity::individual_Ho(dms$gt, genetic_group_var=NULL, maf=0.02, max_missingness = 1)
pca_hybrids$Ho <- ind_ho[match(pca_hybrids$Row.names, names(ind_ho))]

pc1_ho_plot <- ggplot(pca_hybrids,
       aes(x=PC1, y=Ho, colour=morphid2, shape=site))+
  xlab(pcnames[1])+ylab("Ho (MAF 2%)")+
  geom_point(size=1)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+
  geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="HBDSCAN Cluster", shape="Site")+
  theme(legend.key.size = unit(0, 'lines'))+
  guides(colour = guide_legend(title.position = "top"),
         shape = guide_legend(title.position = "top"))+
  scale_colour_manual(values=morphid_colours)

pc1_ho_plot

parent1 <- pca_hybrids$Row.names[which(pca_hybrids$PC1 < -5)]
parent2 <- pca_hybrids$Row.names[which(pca_hybrids$PC1 > 10)]


d3.1_hybrids <- remove.by.list(d3,hybrids_only)
d4_hybrids        <- remove.poor.quality.snps(d3.1_hybrids, min_repro=0.98, max_missing=0.9)%>% remove.fixed.snps()
d5_hybrids        <- sample.one.snp.per.locus.random(d4_hybrids, seed=12345) 

gt_hybrids <- d3$gt[which(d3$sample_names %in% hybrids_only),]


# Calculate allele frequencies with missing data handled
parent1_af <- colSums(gt_hybrids[parent1, ], na.rm = TRUE) / (2 * colSums(!is.na(gt_hybrids[parent1, ])))
parent2_af <- colSums(gt_hybrids[parent2, ], na.rm = TRUE) / (2 * colSums(!is.na(gt_hybrids[parent2, ])))

# Identify loci with high missingness (>90%) in each parent group
parent1_na <- colnames(gt_hybrids[, which((colSums(is.na(gt_hybrids[parent1,])) / length(parent1)) > 0.95)])
parent2_na <- colnames(gt_hybrids[, which((colSums(is.na(gt_hybrids[parent2,])) / length(parent2)) > 0.95)])

# Identify loci fixed for 'a' and 'b' alleles in each parent group
parent1_fixed_a <- setdiff(colnames(gt_hybrids[, which(parent1_af <= 0.05)]), parent1_na)
parent1_fixed_b <- colnames(gt_hybrids[, which(parent1_af >= 0.95)])

parent2_fixed_a <- setdiff(colnames(gt_hybrids[, which(parent2_af <= 0.05)]), parent2_na)
parent2_fixed_b <- colnames(gt_hybrids[, which(parent2_af >= 0.95)])

# Find loci with opposite fixed alleles
keep_loci <- c(parent1_fixed_a[parent1_fixed_a %in% parent2_fixed_b],
               parent1_fixed_b[parent1_fixed_b %in% parent2_fixed_a])

# Subset and remove loci fixed in hybrids
gt_hybrids_fixed <- gt_hybrids[, keep_loci, drop = FALSE]

pca_hybrids <- pca_hybrids[match(rownames(gt_hybrids_fixed), pca_hybrids$Row.names),]
parent_ann <- rowAnnotation(
  PC1 = as.numeric(pca_hybrids$PC1),
  col = list(
    PC1 = colorRamp2(
      c(min(pca_hybrids$PC1, na.rm = TRUE), mean(pca_hybrids$PC1, na.rm = TRUE), max(pca_hybrids$PC1, na.rm = TRUE)),
      c("red", "white", "blue")
    )
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 0)
)

# Plot heatmap
hybrid_gt_heatmap <- Heatmap(as.matrix(gt_hybrids_fixed),
        # col=colorRamp2(c(0, 1, 2), c("#FFEE8C", "forestgreen", "lightblue")),
        col=c("0" = "#FFEE8C", "1" = "forestgreen", "2" = "lightblue"),
        na_col = "white",
        row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = parent_ann,
        heatmap_legend_param = list(
          title = "Genotype",               # Change legend title
          at = c(0, 1, 2),                  # Specify values in the legend
          labels = c("AA", "AB", "BB")      # Custom labels for the legend
        )
        )



combined_plots <- multi_panel_figure(
  width = c(10, 18),   # Adjust these dimensions as needed
  height = c(11,0.01),
  unit = "cm",
  panel_label_type = "upper-roman"
)

# Fill the panels with the respective plots
combined_plots %<>%
  fill_panel(pc1_ho_plot + theme(legend.position = "bottom", legend.direction = 'vertical'), column = 1, row = 1, label = "a") %<>%
  fill_panel(draw(hybrid_gt_heatmap), column = 2, row = 1, label = "b")

ggsave('LantCama/outputs/SupFig_hybrids.pdf', combined_plots, width = 29, height = 12, units = "cm")

