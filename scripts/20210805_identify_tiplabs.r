# Install packages
pacman::p_load("phylogram")
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
torem <- readLines("data/20210601_spliceases_cut_from_dendro.txt") 
newtr <- keep.tip(tr, grep(paste0(torem, collapse = "|"), tr$tip.label))
dend2 <- as.dendrogram.phylo(newtr)

# Try cutting into 20 groups
dendcut20 <- cutree(dend2, k = 20)
