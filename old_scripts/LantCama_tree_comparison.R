# Robinson-Foulds Distance: This method measures the dissimilarity between two trees by counting the number of splits (bipartitions) that are present in one tree but not the other.
# 
# Kendall-Colijn Distance: This method measures the distance between trees by comparing the path lengths between pairs of tips in the trees.

rf_distance1 <- TreeDist::RobinsonFoulds(iqtree, as.phylo(mrbayes))
kc_distance1 <- TreeDist::PathDist(iqtree, as.phylo(mrbayes))

rf_distance1
kc_distance1

rf_distance2 <- TreeDist::RobinsonFoulds(iqtree, upgma)
kc_distance2 <- TreeDist::PathDist(iqtree, upgma)

rf_distance2
kc_distance2


rf_distance3 <- TreeDist::RobinsonFoulds(as.phylo(mrbayes), upgma)
kc_distance3 <- TreeDist::PathDist(as.phylo(mrbayes), upgma)

rf_distance3
kc_distance3

num_leaves <- length(upgma$tip.label)
max_rf_distance <- 2 * (num_leaves - 3)
normalized_rf_distance1 <- rf_distance1 / max_rf_distance
normalized_rf_distance2 <- rf_distance2 / max_rf_distance
normalized_rf_distance3 <- rf_distance3 / max_rf_distance

print(normalized_rf_distance)
# In summary, a "good" RF value depends on your specific context and the size of your trees. Generally, lower values indicate more similar trees, and normalizing the RF distance can help provide a clearer interpretation.


kc_distance1
kc_distance2
kc_distance3
