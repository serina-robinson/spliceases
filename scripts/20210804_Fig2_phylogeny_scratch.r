# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer", "tidytree",
               "Biostrings", "treeio", "ggtree", "scales", "adephylo")

# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
torem <- readLines("data/20210601_spliceases_cut_from_dendro.txt") 
newtr <- keep.tip(tr, grep(paste0(torem, collapse = "|"), tr$tip.label))
midtr <- midpoint.root(newtr)
gt2 <- ggtree(midtr) %>% groupClade(.node=c(1688, 
                                  1829,
                                  1145,
                                  1688,
                                  1288,
                                  986,
                                  1862,
                                  1673,
                                  1644,
                                  1652))  +
  aes(color = group)
gt2

pdf("output/ggtree_grouped.pdf", width = 20, height = 20)
gt3 <- gt2 + 
  aes(color=group, linetype=group)
gt3
dev.off()

gt3
#   ggtree::geom_hilight(node = 1829, fill = "#6b68d4", alpha = .4) + # Frankia CLADE VIII
#   ggtree::geom_hilight(node = 1820, fill = "#55bec5", alpha = .4) + 
#   ggtree::geom_hilight(node = 1145, fill = hex_values[1])+ #alpha = .4) + # KYG clade
#   ggtree::geom_hilight(node = 1688, fill = "forestgreen", alpha = .4) + # 1651 is too far back Cyanobacteria
#   ggtree::geom_hilight(node = 1288, fill = hex_values[6]) +#, alpha = .4) + # GDL Burkholderia 
#   ggtree::geom_hilight(node = 986, fill="maroon", alpha = .4) + # KYG clade 2
#   ggtree::geom_hilight(node = 1862, fill = hex_values[5], alpha = .4) + #Myxobacteria 1
#   ggtree::geom_hilight(node = 1673, fill = hex_values[5], alpha = .4) + #Myxobacteria 2 (Clade V)
#   ggtree::geom_hilight(node = 1644, fill = hex_values[4]) + # NHLP!  +# alpha = .4) + 
#   ggtree::geom_hilight(node = 1652, fill = "firebrick", alpha = .4) + # Pseudoalteromonas Clade II

gtr <- ggtree(tr)
alltips <- adephylo::listTips(tr)
gc <- groupClade(gtr, .node = 1688)

#gc2 <- groupClade(gtr, .node = 1145)

gtree <- ggtree::groupClade(tr, .node=c(1688))

ggtree(gtree, aes(color=group, linetype=group)) 

# dend <- read_csv("data/dendrogram_cuts/20200911_dendrogram_cutting_20grps.csv")
# colnames(dend)[2] <- "dendcut"
# singletons <- as.numeric(names(table(dend$dendcut)[table(dend$dendcut) <= 100]))
# dendr <- dend %>%
#   dplyr::mutate(dendrsplit = case_when(dendcut %in% singletons ~ "singleton",
#                                        TRUE ~ as.character(dendcut)))
# dendr

# Nodes to split



groupInfo <- split(tr$tip.label, dendr$dendrsplit)
groupInfo

phylip <- groupOTU(tr, groupInfo)
gtfake <- ggtree(phylip)

# Use the dendrogram to color
pal <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
pal2 <- brewer.pal(8,"Set1")[c(1:3, 5)]
# pal3 <- c("gray40", pal[c(2)], pal[1],  "#E41A1C", pal[5], pal[6],  "gray70", pal[3], pal[4], "goldenrod4")
pal3 <- c(rep("gray50", 3),  "#E41A1C", rep("gray50", 6))


# ggtree
midtr2 <- midpoint.root(tr)
gtr <- ggtree(midtr2, #color = pal3[as.factor(gtfake$data$group)],
              layout = "fan", open.angle = 0, 
              branch.length = "none", size = 0.25) %>%
              groupClade(.node=c(1688, 
                                       1829,
                                       1145,
                                       1688,
                                       1288,
                                       986,
                                       1862,
                                       1673,
                                       1644,
                                       1652)) +
  aes(color=group)

gtr

metadat <- data.frame(label = gtr$data$label[gtr$data$isTip]) %>%
  dplyr::mutate(query_acc = case_when(grepl("WP_|NP_|YP_", label) ~ 
                                        paste0(word(label, 1, sep = "_"),
                                               "_", word(label, 2, sep = "_")),
                                      TRUE ~ paste0(word(label, 1, sep = "_")))) %>%                                
  dplyr::mutate(id = paste0(word(label, sep = "_", -2), "_", word(label, sep = "_", -1))) %>%
  dplyr::mutate(classification = word(label, sep = "_", -1)) %>%
  dplyr::mutate(phylum = paste0(word(label, sep = "_", 3))) %>%
  dplyr::mutate(class = paste0(word(label, sep = "_", 4))) %>%
  dplyr::mutate(order = paste0(word(label, sep = "_", 5))) %>%
  dplyr::mutate(family = paste0(word(label, sep = "_", 6))) %>%
  dplyr::mutate(genus = paste0(word(label, sep = "_", 7))) %>%
  dplyr::mutate(is_anchor = dplyr::case_when(classification == "PFAM" ~ FALSE,
                                             classification == "characterized" ~ FALSE,
                                             TRUE ~ TRUE)) %>%
  dplyr::mutate(is_interesting = dplyr::case_when(classification == "PFAM" ~ FALSE,
                                                  TRUE ~ TRUE)) %>%
  dplyr::mutate(to_color = dplyr::case_when(classification == "PFAM" ~ "white",
                                            classification == "characterized" ~ "#E31A1C",
                                            TRUE ~ "gray40"))

# Combine metadata with tree
gdf <- gtr %<+% metadat

# Set neighborhood size
# Radical SAM is usually 6
# nbsize <- c(5, 7) 
# nbsize <- c(4, 5, 7, 8)
nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)

# Read in the data for XYG density plotting
gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv")

gbarmerg <- gbardat %>%
  dplyr::distinct(., query, variable, .keep_all =T) %>%
  dplyr::filter(row_id %in% nbsize) %>% # set neighborhood size
  add_count(query_acc, name = "yg_count") %>%
  arrange(genus_species) %>%
  dplyr::mutate(aa_len = nchar(aa)) %>%
  group_by(query) %>%
  dplyr::mutate(yg_density = sum(yg_n_recount)/sum(aa_len))

summary(gbarmerg$yg_density)

mergdat <- metadat %>%
  dplyr::left_join(., gbarmerg, by = "query_acc")

# Now try constraining the neighborhood size
nbdat <- mergdat %>%
  distinct(query_acc, .keep_all = T) %>%
  dplyr::mutate(barcol = ifelse(grepl("characterized",label), "gray40", "#E31A1C"))

nbdat$yg_count[is.na(nbdat$yg_count)] <- 0
nbdat$yg_density[is.na(nbdat$yg_density)] <- 0
sort(table(nbdat$yg_density))

# Set color palette
pal <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
pal2 <- brewer.pal(8,"Set1")[c(1:3, 5)]
pal2
pal3 <- c("gray40", pal[2],"red", "#8DA0CB", #"#377EB8",
          pal[c(5, 4)], "gray40", pal[6],  "navyblue","goldenrod4")

pdf(paste0("output/20210803_yg_density_", length(nbsize), "genes.pdf"), width = 10, height = 10)
gpl <- gdf + 
  ggtreeExtra::geom_fruit(
    data = nbdat,
    geom = geom_bar, 
    mapping = aes(y = label, x = yg_density, fill = barcol),
    axis.params=list(
      axis   = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5
    ),
    grid.params=list(),
    orientation = "y",
    pwidth = 0.5,
    offset = 0.01,
    stat = "identity") +
  scale_fill_manual(values = pal3) +
  theme(legend.position = "none")
gpl
dev.off()


