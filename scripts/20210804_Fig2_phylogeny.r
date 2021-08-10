# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer",
               "Biostrings", "treeio", "ggtree", "scales")

# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")

dend <- read_csv("data/dendrogram_cuts/20200911_dendrogram_cutting_20grps.csv")
colnames(dend)[2] <- "dendcut"
singletons <- as.numeric(names(table(dend$dendcut)[table(dend$dendcut) <= 100]))
dendr <- dend %>%
  dplyr::mutate(dendrsplit = case_when(dendcut %in% singletons ~ "singleton",
                                       TRUE ~ as.character(dendcut)))

groupInfo <- split(tr$tip.label, dendr$dendrsplit)
phylip <- groupOTU(tr, groupInfo)
gtfake <- ggtree(phylip)

# Use the dendrogram to color
pal <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
pal2 <- brewer.pal(8,"Set1")[c(1:3, 5)]
# pal3 <- c("gray40", pal[c(2)], pal[1],  "#E41A1C", pal[5], pal[6],  "gray70", pal[3], pal[4], "goldenrod4")
pal3 <- c(rep("gray50", 3),  "#E31A1C", rep("gray50", 6))


# ggtree
gtr <- ggtree(tr, color = pal3[as.factor(gtfake$data$group)],
              layout = "fan", open.angle = 0, 
              branch.length = "none", size = 0.25)

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
# bsize <- c(4, 5, 7, 8)
nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)

# Read in the data for XYG density plotting
gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv")

gbarmerg <- gbardat %>%
  #dplyr::distinct(., query, variable, .keep_all =T) %>%
  dplyr::mutate(uniq_id = paste0(protein_acc, "_", variable)) %>%
  dplyr::filter(!duplicated(uniq_id)) %>%
  dplyr::filter(row_id %in% nbsize) %>% # set neighborhood size
  add_count(query_acc, name = "yg_count") %>%
  arrange(genus_species) %>%
  dplyr::mutate(aa_len = nchar(aa)) %>%
  group_by(query) %>%
  dplyr::mutate(yg_density = sum(yg_n_recount)/sum(aa_len)) %>%
  ungroup()

summary(gbarmerg$yg_density)

mergdat <- metadat %>%
  dplyr::left_join(., gbarmerg, by = "query_acc")

# Now try constraining the neighborhood size
splics <- readLines("data/20210601_spliceases_cut_from_dendro.txt") 
label %in% splics

nbdat <- mergdat %>%
  arrange(desc(yg_density)) %>%
  distinct(query_acc, .keep_all = T) %>%
  dplyr::mutate(barcol = ifelse(label %in% splics, "#E31A1C",  "gray40"))

nbdat$yg_count[is.na(nbdat$yg_count)] <- 0
nbdat$yg_density[is.na(nbdat$yg_density)] <- 0
sort(table(nbdat$yg_density))

# Set color palette
pal <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
pal2 <- brewer.pal(8,"Set1")[c(1:3, 5)]
pal2
pal3 <- c("#E31A1C", "gray40", pal[2],"red", "#8DA0CB", #"#377EB8",
          pal[c(5, 4)], "gray40", pal[6],  "navyblue","goldenrod4")
table(nbdat$barcol)


pdf(paste0("output/20210803_yg_density_", length(nbsize), "_arranged_genes.pdf"), width = 10, height = 10)
gpl <- gdf + 
  ggtreeExtra::geom_fruit(
    data = nbdat,
    geom = geom_bar, 
    mapping = aes(y = label, x = yg_density, fill = barcol, color = barcol), 
    axis.params=list(
      line.color = "black",
      axis   = "x",
      text.size  = 4,
      hjust      = 1,
      vjust      = 0.5
    ),
    grid.params=list(color = "gray70"),
    orientation = "y",
    pwidth = 0.5,
    offset = 0.01,
    stat = "identity") +
  scale_fill_manual(values = pal3) +
  scale_color_manual(values = pal3) +
  theme(legend.position = "none")
rotate_tree(gpl, 90)
dev.off()

table(nbdat$barcol)
