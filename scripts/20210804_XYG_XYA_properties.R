# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer",
               "Biostrings", "treeio", "ggtree", "scales")

# Read in the data for XYG density plotting
gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv")

# Keep only the putative spliceases
splics <- readLines("data/20210601_spliceases_cut_from_dendro.txt")
seqselect <- tr$tip.label[tr$tip.label %in% splics]
querselect <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
                          paste0(word(seqselect, 1, sep = "_"),
                                 "_", word(seqselect, 2, sep = "_")),
                        TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# Set neighborhood size
# Radical SAM is usually 6
# nbsize <- c(5, 7) 
# bsize <- c(4, 5, 7, 8)
nbsize <- c(3, 4, 5, 7, 8, 9)
# nbsize <- c(2, 3, 4, 5, 7, 8, 9, 10) 
# nbsize  <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)

# Filter only for spliceases
gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv") %>%
    dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) 
  
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