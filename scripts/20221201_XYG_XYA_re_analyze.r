# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer", "patchwork", "microshades",
               "Biostrings", "treeio", "ggtree", "scales", "ggpubr")

# Keep only the putative spliceases
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
splics <- readLines("data/20210601_spliceases_cut_from_dendro.txt")
seqselect <- tr$tip.label[tr$tip.label %in% splics]
querselect <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
                          paste0(word(seqselect, 1, sep = "_"),
                                 "_", word(seqselect, 2, sep = "_")),
                        TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# Set neighborhood size of plus/minus 3 ORFs
nbsize <- c(3, 4, 5, 7, 8, 9)

# Filter only for spliceases
gbardat <- read_csv("data/XYG_count_data/unfiltered_combined_YG_YA_count_data.csv") %>%
    dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) %>%
    dplyr::filter(row_id %in% nbsize) %>%
    dplyr::mutate(aa_class = case_when(value %in% c("G", "A", "V", "L", "I",
                                                    "M", "F", "Y", "P", "C", "W") ~ "nonpolar",
                                       value %in% c("S", "T", "Y", "N", "Q") ~
                                                      "polar",
                                       value %in% c("R", "K", "H", "D", "E") ~ "charged")) 
yg_dat <- gbardat %>%
  dplyr::filter(grepl("yg_", variable)) 
relevel_sort <- names(sort(table(yg_dat$value), decreasing = T))
yg_dat$value <- factor(yg_dat$value, levels = relevel_sort) # to order bars in bar plot

ya_dat <- gbardat %>%
  dplyr::filter(grepl("ya_", variable))  
ya_dat$value <- factor(ya_dat$value, levels = relevel_sort) # to order bars in bar plot

ya_dat_check <- ya_dat %>%
  dplyr::mutate(aa_len = width(aa)) %>%
  dplyr::select(value, row_id, query, genus_species, nucleotide_acc, protein_acc, name1,
                pfam_id1, description1, variable, aa, aa_len)
write_csv(ya_dat_check, "output/ya_plus_mins_3_ORFs_check_aa.csv")

yg_dat_check <- yg_dat %>%
  dplyr::mutate(aa_len = width(aa)) %>%
  dplyr::select(value, row_id, query, genus_species, nucleotide_acc, protein_acc, name1,
                pfam_id1, description1, variable, aa, aa_len)
write_csv(yg_dat_check, "output/yg_plus_mins_3_ORFs_check_aa.csv")

# Probe the number of unique XYA and XYG hits
overlap <- length(intersect(unique(yg_dat$protein_acc), unique(ya_dat$protein_acc)))
xyg_unique <- length(unique(yg_dat$protein_acc)) 
xya_unique <- length(unique(ya_dat$protein_acc))

# Calculate the substrate properties
yg_max_3 <- yg_dat %>%
  dplyr::mutate(aa_len = width(aa)) %>%
  dplyr::filter(row_id %in% c(3,4,5,7,8,9)) %>%
  dplyr::filter(!duplicated(protein_acc)) %>%
  dplyr::arrange(aa_len) %>%
  dplyr::filter(aa_len >= 600) %>%
  dplyr::select(row_id, query, genus_species, nucleotide_acc, protein_acc, name1,
                pfam_id1, description1, aa, aa_len)
write_csv(yg_max_3, "output/longest_substrates_xyg_plus_minus_3_ORFs.csv")

yg_max_1 <- yg_dat %>%
  dplyr::mutate(aa_len = width(aa)) %>%
  dplyr::filter(row_id %in% c(5,7)) %>%
  dplyr::filter(!duplicated(protein_acc)) %>%
  dplyr::arrange(aa_len) %>%
  dplyr::slice(1:20, (nrow(.)-19):nrow(.)) %>%
  dplyr::select(row_id, query, genus_species, nucleotide_acc, protein_acc, name1,
                pfam_id1, description1, aa, aa_len)
write_csv(yg_max_1, "output/xyg_plus_minus_1_ORFs.csv")


# Set color palette
hex_values <-c(microshades_palette("micro_purple", 4, lightest = TRUE)[4], 
               microshades_palette("micro_blue", 3, lightest = TRUE)[3],
               microshades_palette("micro_orange", 3, lightest = TRUE)[3])

ygbar_6 <- ggplot(yg_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0))
ygbar_6 <- ygbar_6 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

yabar_6 <- ggplot(ya_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0))
yabar_6 <- yabar_6 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) 

# Set neighborhood size plus/minus 1 ORF
 nbsize <- c(5, 7) 

 gbardat <- read_csv("data/XYG_count_data/unfiltered_combined_YG_YA_count_data.csv") %>%
   dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) %>%
   dplyr::filter(row_id %in% nbsize) %>%
   dplyr::mutate(aa_class = case_when(value %in% c("G", "A", "V", "L", "I",
                                                   "M", "F", "Y", "P", "C", "W") ~ "nonpolar",
                                      value %in% c("S", "T", "Y", "N", "Q") ~
                                        "polar",
                                      value %in% c("R", "K", "H", "D", "E") ~ "charged")) 
 

 yg_dat <- gbardat %>%
   dplyr::filter(grepl("yg_", variable)) 
 relevel_sort <- names(sort(table(yg_dat$value), decreasing = T))
 yg_dat$value <- factor(yg_dat$value, levels = relevel_sort)
 
 ya_dat <- gbardat %>%
   dplyr::filter(grepl("ya_", variable)) 
 ya_dat$value <- factor(ya_dat$value, levels = relevel_sort)
 
# Set color palette
hex_values <-c(microshades_palette("micro_purple", 4, lightest = TRUE)[4], 
               microshades_palette("micro_blue", 3, lightest = TRUE)[3],
               microshades_palette("micro_orange", 3, lightest = TRUE)[3])

# Bar plot
ygbar_2 <- ggplot(yg_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0))
ygbar_2 <- ygbar_2 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank())

yabar_2 <- ggplot(ya_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0, 0))

yabar_2 <- yabar_2 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank())
yabar_2

layout <- '
AABB
CCDD'


# Print out everything in one PDF
pdf("output/unfiltered_supplemental_YA_YG.pdf", width = 8, height = 8)
wrap_plots(
  A = ygbar_2,
  B= ygbar_6,
  C = yabar_2,
  D = yabar_6,
  design = layout)
dev.off()
