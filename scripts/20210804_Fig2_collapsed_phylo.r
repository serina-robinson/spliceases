# Install packages
# remotes::install_github("KarstensLab/microshades")
pacman::p_load("ggtree", "stringr", "ggrepel", "dplyr", "scales", "phytools",
               "ggtreeExtra", "ggnewscale", "readr", "RColorBrewer", "microshades")

# Set color-blind friendly color palette
hex_values <-c(microshades_palette("micro_cvd_purple", 1, lightest = TRUE), 
               microshades_palette("micro_cvd_blue", 2, lightest = TRUE)[2],
               microshades_palette("micro_cvd_green", 2, lightest = TRUE)[2],
               microshades_palette("micro_cvd_turquoise", 2, lightest = TRUE)[2],
               microshades_palette("micro_brown", 2, lightest = TRUE)[2],
               microshades_palette("micro_orange", 1, lightest = TRUE),
               microshades_palette("micro_cvd_gray", 3, lightest = TRUE)[3])
hex_values

pal2 <- c("#6b68d4",
          "#55bec5",
          "#EFB6D6",
          "forestgreen",
          "#feeda0",
          "#CAA995",
          "#48C9B0",
          "#B7B7B7")


# Read in tree and trim
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
torem <- readLines("data/20210601_spliceases_cut_from_dendro.txt") 
newtr <- keep.tip(tr, grep(paste0(torem, collapse = "|"), tr$tip.label))
length(torem)
write.tree(newtr, file = "data/seqs_for_phylogeny/966_seqs_splicease_phylo.nwk")
midtr <- midpoint.root(newtr)
write.tree(midtr, file = "data/seqs_for_phylogeny/966_seqs_splicease_phylo_rooted_at_midpoint.nwk")

# Clades to collapse

# Clade 968 should be dropped from the tree # Or an outgroup
# 1651 is Cyanobacteria + Pseudoalteromonas
# 1355 is GDL Burkholderia
# 1830/1831 is Frankia
# 1096 and 986 are the KYG clade
gdf <- ggtree(midtr, layout = "rectangular", color = "#E41A1C", size = 1)

# New colors
pdf(paste0("output/trees/20210809_spliceases_only_full_labels_with_archaea.pdf"), width = 21.59, height = 27.94)
gdf + 
  geom_tiplab(size = 0.8) +
  geom_nodelab(aes(label = node), size = 0.8) +
  ggtree::geom_hilight(node = 1829, fill = "#6b68d4", alpha = .4) + # Frankia CLADE VIII
  ggtree::geom_hilight(node = 1820, fill = "#55bec5", alpha = .4) + 
  ggtree::geom_hilight(node = 1096, fill = "dodgerblue")+ #alpha = .4) + # KYG clade
  ggtree::geom_hilight(node = 1688, fill = "forestgreen", alpha = .4) + # 1651 is too far back Cyanobacteria
  ggtree::geom_hilight(node = 1288, fill = hex_values[6], alpha = .4) + # GDL Burkholderia 
  ggtree::geom_hilight(node = 986, fill="maroon", alpha = .4) + # KYG clade 2
  ggtree::geom_hilight(node = 1862, fill = hex_values[5], alpha = .4) + #Myxobacteria 1
  ggtree::geom_hilight(node = 1673, fill = hex_values[5], alpha = .4) + #Myxobacteria 2 (Clade V)
  ggtree::geom_hilight(node = 1644, fill = hex_values[4]) + # NHLP!  +# alpha = .4) + 
  ggtree::geom_hilight(node = 1652, fill = hex_values[1], alpha = .4) + # Pseudoalteromonas Clade II
  ggtree::geom_hilight(node = 1232, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1110, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1635, fill = hex_values[7], alpha = .4) + # uncharacterized Bradymonadales
  ggtree::geom_hilight(node = 1844, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1919, fill = hex_values[7], alpha = .4) + # uncharacterized sister clade to Myxos
  ggtree::geom_hilight(node = 986, fill = hex_values[7], alpha = .4) + # uncharacterized Bradymonadales
  ggtree::geom_hilight(node = 968, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1919, fill = hex_values[7], alpha = .4) + # uncharacterized sister clade to Myxos
  theme(legend.position = "none") +
  xlim_tree(5)
dev.off()


pdf(paste0("output/20210809_spliceases_only_full_labels_with_archaea.pdf"), width = 21.59, height = 27.94)
gdf + 
  geom_tiplab(size = 0.8) +
  geom_nodelab(aes(label = node), size = 0.8) +
  ggtree::geom_hilight(node = 1829, fill = "#6b68d4", alpha = .4) + # Frankia CLADE VIII
  ggtree::geom_hilight(node = 1820, fill = "#55bec5", alpha = .4) + 
  ggtree::geom_hilight(node = 1108, fill = hex_values[1]) + #alpha = .4) + # KYG clade #1145 was proper
  ggtree::geom_hilight(node = 1688, fill = "forestgreen", alpha = .4) + # 1651 is too far back Cyanobacteria
  ggtree::geom_hilight(node = 1288, fill = hex_values[6]) +#, alpha = .4) + # GDL Burkholderia 
  ggtree::geom_hilight(node = 986, fill="maroon", alpha = .4) + # KYG clade 2
  ggtree::geom_hilight(node = 1862, fill = hex_values[5], alpha = .4) + #Myxobacteria 1
  ggtree::geom_hilight(node = 1673, fill = hex_values[5], alpha = .4) + #Myxobacteria 2 (Clade V)
  ggtree::geom_hilight(node = 1644, fill = hex_values[4]) + # NHLP!  +# alpha = .4) + 
  ggtree::geom_hilight(node = 1652, fill = "firebrick", alpha = .4) + # Pseudoalteromonas Clade II
  ggtree::geom_hilight(node = 1635, fill = hex_values[7], alpha = .4) + # uncharacterized Bradymonadales
  ggtree::geom_hilight(node = 1844, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1919, fill = hex_values[7], alpha = .4) + # uncharacterized sister clade to Myxos
  ggtree::geom_hilight(node = 986, fill = hex_values[7], alpha = .4) + # uncharacterized Bradymonadales
  ggtree::geom_hilight(node = 968, fill = hex_values[7], alpha = .4) + # uncharacterized
  ggtree::geom_hilight(node = 1097, fill = hex_values[7], alpha = .4) + 
  # ggtree::geom_hilight(node = 1919, fill = hex_values[7], alpha = .4) + 
 # ggtree::geom_hilight(node = 1205, fill = hex_values[7], alpha = .4) + # uncharacterized sister clade to Myxos
  theme(legend.position = "none") +
  xlim_tree(5)
dev.off()


pdf(paste0("output/20210809_spliceases_collapsed.pdf"), width = 21.59, height = 27.94)
gdf %>%
  ggtree::collapse(., 1829, 'min', fill = "#6b68d4", alpha = .4) %>%
  ggtree::collapse(., 1820, 'min', fill = "#55bec5", alpha = .4) %>%
  ggtree::collapse(., 1136,  'min', fill = "dodgerblue", alpha = .4)  %>% #alpha = .4) + # KYG clade #1145 was proper
  ggtree::collapse(., 1688,  'min', fill = "forestgreen", alpha = .3) %>% # 1651 is too far back Cyanobacteria
  ggtree::collapse(., 1288,  'min', fill = hex_values[6], alpha = .4) %>% # GDL Burkholderia 
  ggtree::collapse(., 986, 'min', fill = "orange", alpha = .5) %>%  # KYG clade 2
  ggtree::collapse(., 1862, 'min', fill = hex_values[5], alpha = .4) %>% #Myxobacteria 1
  ggtree::collapse(., 1673, 'min', fill = hex_values[5], alpha = .4) %>% #Myxobacteria 2 (Clade V)
  ggtree::collapse(., 1644, 'min', fill = hex_values[4]) %>% # NHLP!  +# alpha = .4) + 
  ggtree::collapse(., 1652, 'min', fill = hex_values[1], alpha = .4) %>% # Pseudoalteromonas Clade II
  ggtree::collapse(., 1635, 'min', fill = hex_values[7], alpha = .4) %>% # uncharacterized Bradymonadales
  ggtree::collapse(., 1844, 'min', fill = hex_values[5], alpha = .4) %>% # uncharacterized
  ggtree::collapse(., 1919, 'min', fill = hex_values[7], alpha = .4) %>% # uncharacterized sister clade to Myxos
  ggtree::collapse(., 968, 'min', fill = hex_values[7], alpha = .4) %>% # uncharacterized
  ggtree::collapse(., 1097, 'min', fill = hex_values[7], alpha = .4) %>%
  ggtree::collapse(., 1686, 'min', fill = hex_values[7], alpha = .4) %>% # uncharacterized
  ggtree::collapse(., 1232, 'min', fill = hex_values[7], alpha = .4) %>% # uncharacterized
  ggtree::collapse(., 1110, 'min', fill = hex_values[7], alpha = .4) + # uncharacterized
  # ggtree::geom_hilight(node = 1919, fill = hex_values[7], alpha = .4) + 
  # ggtree::geom_hilight(node = 1205, fill = hex_values[7], alpha = .4) + # uncharacterized sister clade to Myxos
  theme(legend.position = "none") +
  xlim_tree(1)
dev.off()



