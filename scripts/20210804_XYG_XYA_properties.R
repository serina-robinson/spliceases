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

# Set neighborhood size
# Radical SAM is usually 6
nbsize <- c(5, 7) 
# nbsize <- c(4, 5, 7, 8)
# nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)



# Filter only for spliceases
gbardat <- read_csv("data/XYG_count_data/combined_YG_YA_count_data.csv") %>%
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
ya_dat

# Set color palette
hex_values <-c(microshades_palette("micro_purple", 4, lightest = TRUE)[4], 
               microshades_palette("micro_blue", 3, lightest = TRUE)[3],
               microshades_palette("micro_orange", 3, lightest = TRUE)[3])
#pal2 <- c("#E0DAEA", hex_values[2], hex_values[3])


# Bar plot
ygbar_2 <- ggplot(yg_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values)
ygbar_2 <- ygbar_2 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank())# +
  #ggtitle('A')


yabar_2 <- ggplot(ya_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values)
yabar_2 <- yabar_2 +
  #ggtitle('C') +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank())


#### NB size of 6
# Set neighborhood size
# Radical SAM is usually 6
# nbsize <- c(5, 7) 
# nbsize <- c(4, 5, 7, 8)
nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)



# Filter only for spliceases
gbardat <- read_csv("data/XYG_count_data/combined_YG_YA_count_data.csv") %>%
  dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) %>%
  dplyr::filter(row_id %in% nbsize) %>%
  dplyr::mutate(aa_class = case_when(value %in% c("G", "A", "V", "L", "I",
                                                  "M", "F", "Y", "P", "C") ~ "nonpolar",
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
ya_dat

# Set color palette
hex_values <-c(microshades_palette("micro_purple", 4, lightest = TRUE)[4], 
               microshades_palette("micro_blue", 3, lightest = TRUE)[3],
               microshades_palette("micro_orange", 3, lightest = TRUE)[3])
ygbar_6 <- ggplot(yg_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values)
ygbar_6 <- ygbar_6 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank())


yabar_6 <- ggplot(ya_dat, aes(x = value)) +
  geom_bar(aes(fill = aa_class), position = "dodge", stat = "count") +
  theme_pubr() +
  scale_fill_manual(values = hex_values)
yabar_6 <- yabar_6 +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank())


layout <- '
AABB
CCDD'
# Print out everything in one PDF
pdf("output/supplemental_YA_YG.pdf", width = 8.5, height = 4)
wrap_plots(
  A = ygbar_2,
  B= ygbar_6,
  C = yabar_2,
  D = yabar_6,
  design = layout)
dev.off()
