# Install packages
pacman::p_load("ggtree", "stringr", "ggrepel", "dplyr", "scales", "phytools", "Biostrings",
               "ggtreeExtra", "ggnewscale", "readr", "RColorBrewer", "DECIPHER")

# Read in the consensus tree
# midtr <- read.tree("data/seqs_for_phylogeny/2814_fasttree_gt_10_rooted_at_midpoint.nwk")

#tr <- read.tree("data/seqs_for_phylogeny/2814_fasttree_gt_10.nwk")
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
torem <- readLines("output/20210601_spliceases_cut_from_dendro.txt") 
length(torem) # 966 spliceases
newtr <- keep.tip(tr, grep(paste0(torem, collapse = "|"), tr$tip.label))
midtr <- midpoint.root(newtr)

gtr <- ggtree(midtr,
              layout = "rectangular",
              size = 0.25) 

phys <- readLines("output/phylums_to_color.txt")
phys

metadat <- data.frame(label = gtr$data$label[gtr$data$isTip]) %>%
  # dplyr::mutate(grp = gtr$data$group[gtr$data$isTip]) %>%
  dplyr::mutate(query_acc = case_when(grepl("WP_|NP_|YP_", label) ~ 
                                        paste0(word(label, 1, sep = "_"),
                                               "_", word(label, 2, sep = "_")),
                                      TRUE ~ paste0(word(label, 1, sep = "_")))) %>%                                
  dplyr::mutate(id = paste0(word(label, sep = "_", -2), "_", word(label, sep = "_", -1))) %>%
  dplyr::mutate(classification = word(label, sep = "_", -1)) %>%
  dplyr::mutate(phylum = paste0(word(label, sep = "_", -6))) %>%
  dplyr::mutate(class = paste0(word(label, sep = "_", -5))) %>%
  dplyr::mutate(order = paste0(word(label, sep = "_", -4))) %>%
  dplyr::mutate(family = paste0(word(label, sep = "_", -3))) %>%
  dplyr::mutate(genus = paste0(word(label, sep = "_", -2))) %>%
  dplyr::mutate(color = case_when(
    #phylum == "Unknown" ~ "Other",
    phylum %in% c("III", "Ciliophora", "No") ~ "Other",
    TRUE ~ phylum)) %>%
  dplyr::mutate(is_anchor = dplyr::case_when(classification == "PFAM" ~ FALSE,
                                             classification == "characterized" ~ FALSE,
                                             classification == "BLAST" ~ FALSE,
                                             TRUE ~ TRUE)) %>%
  dplyr::mutate(is_interesting = dplyr::case_when(classification == "PFAM" ~ FALSE,
                                                  classification == "BLAST" ~ FALSE,
                                                  TRUE ~ TRUE)) %>%
  dplyr::mutate(to_color = dplyr::case_when(classification == "PFAM" ~ "gray",
                                            classification == "BLAST" ~ "red",
                                            classification == "characterized" ~ "red",
                                            TRUE ~ "black")) 

# Calculate genus-level statistics
write_csv(data.frame(sort(table(metadat$genus), decreasing = T)), "data/splicease_genera_only.csv")
sort(table(metadat$class), decreasing = T) #Myxococcales
metadat$label[grep("Myxococcales", metadat$label)]
length(grep("Myxococcales", metadat$label))/length(metadat$label) * 100


# Calculate the fraction that are GDL...
gdl <- readLines("output/20210601_GDLs_accs.txt") 
gdl_pulled <- metadat$label[grep(paste0(gdl, collapse = "|"), metadat$label)]
length(gdl_pulled)/length(metadat$label) * 100

# Combine metadata with tree
gtfake <- gtr %<+% metadat
length(gtfake$data$label)
gtr <- ggtree(midtr, #color = pal3[as.factor(gtfake$data$phylum)],
              color = "gray40",
              # color = ifelse((gtfake$data$group) %in% c("49"), "#E41A1C", "gray50"),
              # layout = "fan", open.angle = 0,  branch.length = "none",
              layout = "rectangular", branch.length = "none",
              size = 0.25) 

# Geom_bar_data
gbardat1 <- read_csv("data/motif_data/YG_characterized_barplot_data.csv")
qs1 <- unique(gbardat1$query)

gbardat2 <- read_csv("data/motif_data/YG_1886_nbs_barplot_data.csv")
qs2 <- unique(gbardat2$query)

gbardat3 <- read_csv("data/motif_data/YG_spasm_key_nbs_barplot.csv")
qs3 <- unique(gbardat3$query)

gbardat4 <- read_csv("data/motif_data/YG_blast_key_nbs_barplot.csv")
qs4 <- unique(gbardat4$query)

# Remove duplicates between the 3 datasets
dups1 <- intersect(x = qs1, y = qs2) # none
dups2 <- intersect(x = qs2, y = qs3) # 1
dups3 <- intersect(x = qs1, y = qs3) # lots
dups4 <- intersect(x = qs4, y = qs1)
dups5 <- intersect(x = qs4, y = qs2)
dups6 <- intersect(x = qs4, y = qs3)
dedup <- unique(c(dups1, dups2, dups3, dups4, dups5, dups6))
dedup

# Set neighborhood size
# Radical SAM is usually 6
nbsize <- c(5, 7)  # one gene on either side
# nbsize <- c(4, 5, 7, 8)
# nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)

# Now try constraining the neighborhood size and/or the identity of 'X' in 'XYG'
# x <- c("G","A", "V", "L", "S", "P", "Q", "I", "M")
# "M", "I", "K")
# x <- "G"
# x <- "A"
# x <- "V"
# x <- "L"
# x <- "S"
# x <- "P"
# x <- "Q"
# x <- "I"
# x <- "M"
# x <- "K"

gbarmerg <- gbardat4 %>%
  bind_rows(gbardat1, gbardat2, gbardat3) %>%
  # dplyr::mutate(uniq_id = paste0(query_acc, "_", variable)) %>%
  dplyr::mutate(uniq_id = paste0(protein_acc, "_", variable)) %>%
  dplyr::filter(!duplicated(uniq_id)) %>%
  dplyr::filter(row_id %in% nbsize) %>% # set neighborhood size
  add_count(query_acc, name = "yg_count") %>%
  arrange(genus_species) %>%
  dplyr::mutate(aa_len = nchar(aa)) %>%
  group_by(query) %>%
  dplyr::filter(!duplicated(protein_acc)) %>%
  dplyr::mutate(yg_density = (yg_count/sum(aa_len)) * 100)

mergdat <- metadat %>%
  dplyr::left_join(., gbarmerg, by = "query_acc")

# nadat <- mergdat[is.na(mergdat$n),] # TODO must deal with this!

nbdat <- mergdat %>%
  distinct(query_acc, .keep_all = T) 

nbdat$yg_count[is.na(nbdat$yg_count)] <- 0
nbdat$yg_density[is.na(nbdat$yg_density)] <- 0
sort(table(nbdat$yg_count))
as.factor(gdf$data$color)


gdf <- ggtree(midtr) %<+% metadat

gdf <- ggtree(midtr, #color = pal3[as.factor(gdf$data$color)],
              color = "gray80",
              # layout = "fan", open.angle = 0,  branch.length = "none",
              layout = "rectangular", branch.length = "none",
              size = 0.5) %<+% metadat
gdf$data$color[is.na(gdf$data$color)] <- "gray40"
gdf$data$color
#gdf <- ggtree(midtr, color = pal3[as.factor(gdf$data$color)],
#              layout = "fan", open.angle = 0,  branch.length = "none",
#              size = 0.5) %<+% metadat


pdf(paste0("output/trees/20210601_spliceases_only_no_labels_with_archaea_labels.pdf"), width = 10, height = 10)
gdf +
  #aes(color = pal3[as.factor(gdf$data$color)]) +
  # %
  #  aes(color = I(phylum)) +
  geom_tippoint(size = ifelse(gdf$data$is_interesting[gdf$data$isTip], 1, NA), 
                color = gdf$data$to_color[gdf$data$isTip]) +
  # ggtreeExtra::geom_fruit(
  #   data = nbdat,
  #   geom = geom_bar, 
  #   mapping = aes(y = label, x = yg_density),
  #   fill = "gray40",
  #   axis.params=list(
  #     axis   = "x",
  #     text.size  = 1.8,
  #     hjust      = 1,
  #     vjust      = 0.5
  #   ),
#   grid.params=list(),
#   orientation = "y",
#   pwidth = 0.5,
#   offset = 0.01,
#   stat = "identity") +
scale_fill_manual(values = pal3,
                  na.value = "gray40") +
  theme(legend.position = "none") #+
# geom_label_repel(label = ifelse(gdf$data$label %in% to_label, gdf$data$id, ""),
#                  size = 3, force = 2,
#                  segment.alpha = 1, fill = "white",
#                  segment.color = "gray60",
#                  # box.padding = 1,
#                  nudge_x = 5, nudge_y = 5)
# gpl
# rotate_tree(gpl, 220)
dev.off()

gdf$data$node

pdf(paste0("output/trees/20210601_spliceases_only_full_labels_with_archaea_updated.pdf"), width = 25, height = 20)
gdf + 
  geom_tippoint(size = ifelse(gdf$data$is_interesting[gdf$data$isTip], 1, NA),
                color = gdf$data$to_color[gdf$data$isTip]) +
  #geom_tiplab2(ifelse(gdf$data$to_color[gdf$data$isTip] == "red", gdf$data$label[gdf$data$isTip], NA)) +
  geom_tiplab(size = 0.5) +
  geom_nodelab(aes(label = node)) +
  #xlim(25, NA) +
  # geom_label(aes(label = node)) +
  # ggtreeExtra::geom_fruit(
  #   data = nbdat,
  #   geom = geom_bar, 
  #   mapping = aes(y = label, x = yg_density),
  #   fill = "gray40",
  #   axis.params=list(
  #     axis   = "x",
  #     text.size  = 1.8,
  #     hjust      = 1,
#     vjust      = 0.5
#   ),
#   grid.params=list(),
#   orientation = "y",
#   pwidth = 0.5,
#   offset = 0.01,
#   stat = "identity") +
# scale_fill_manual(values = pal4) +
theme(legend.position = "none") #+
# geom_label_repel(label = ifelse(gdf$data$label %in% to_label, gdf$data$id, ""),
#                  size = 3, force = 2,
#                  segment.alpha = 1, fill = "white",
#                  segment.color = "gray60",
#                  # box.padding = 1,
#                  nudge_x = 5, nudge_y = 5)
# rotate_tree(gpl, 220)
dev.off()


pdf(paste0("output/trees/20210601_spliceases_with_archaea_updated.pdf"), width = 10, height = 10)
gdf + 
  geom_tippoint(size = ifelse(gdf$data$is_interesting[gdf$data$isTip], 1, NA), 
                color = gdf$data$to_color[gdf$data$isTip]) +
  # ggtreeExtra::geom_fruit(
  #   data = nbdat,
  #   geom = geom_bar, 
  #   mapping = aes(y = label, x = yg_density),
  #   fill = "gray40",
  #   axis.params=list(
  #     axis   = "x",
  #     text.size  = 1.8,
  #     hjust      = 1,
  #     vjust      = 0.5
  #   ),
#   grid.params=list(),
#   orientation = "y",
#   pwidth = 0.5,
#   offset = 0.01,
#   stat = "identity") +
# scale_fill_manual(values = pal4) +
theme(legend.position = "none") +
  geom_label_repel(label = ifelse(gdf$data$is_interesting, gdf$data$id, ""),
                   size = 2, force = 2,
                   segment.alpha = 1, fill = "white",
                   segment.color = "gray60",
                   # box.padding = 1,
                   nudge_x = 5, nudge_y = 5)
# rotate_tree(gpl, 220)
dev.off()

txt <- levels(as.factor(gdf$data$color))[!levels(as.factor(gdf$data$color)) %in% c("gray40", "Other")]
pals <- pal3[pal3!="gray40"][1:14]

pdf(paste0("output/phylum_legend_20210601_with_archaea_updated.pdf"),width=5,height=6)
plot.new()
legend("topleft", legend=txt, fill=pals,
       border=FALSE, bty="n", title = "Phylum")
dev.off()
