# Install packages
pacman::p_load("palmerpenguins", "ggstatsplot",
               "ggbeeswarm", "ggpubr", "rsample")

# Read in the data
dat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv")
splic_keep <- readLines("data/20210601_spliceases_cut_from_dendro.txt")
query_rem <- case_when(grepl("WP_|NP_|YP_", splic_keep) ~ 
                       paste0(word(splic_keep, 1, sep = "_"),
                              "_", word(splic_keep, 2, sep = "_")),
                     TRUE ~ paste0(word(splic_keep, 1, sep = "_")))   

# Read in the groups
grps <- read_csv("data/20210805_grouped_tip_labels.csv") %>%
  dplyr::mutate(query_acc = case_when(grepl("WP_|NP_|YP_", label) ~ 
              paste0(word(label, 1, sep = "_"),
                     "_", word(label, 2, sep = "_")),
            TRUE ~ paste0(word(label, 1, sep = "_"))))  

# Join data 
datjoin <- dat %>%
  left_join(., grps, by = "query_acc")

# Compare spliceases to non-spliceases
subdat <- datjoin %>%
  dplyr::mutate(is_splicease = ifelse(query_acc %in% query_rem, "Splicease", "Other rSAM-SPASM")) %>%
  dplyr::mutate(aa_len = nchar(aa)) %>%
  dplyr::mutate(yg_density = (yg_n_recount/sum(aa_len))) %>%
  dplyr::mutate(uniq_id = paste0(protein_acc, "_", variable)) %>%
  dplyr::filter(!duplicated(uniq_id)) %>%
  dplyr::filter(yg_n_recount != 0)
datjoin$group[is.na(datjoin$group)] <- 11


# plt2 <- ggstatsplot::ggbetweenstats(
#   data = subdat,
#   x = is_splicease,
#   y = yg_n_recount
# )
# plt2

# Test a random sub-sample
set.seed(1234)
dat_split <- initial_split(subdat, prop = 1/10, strata = "is_splicease")
dat_train <- training(dat_split)
dat_train

datjoin$group <- as.factor(datjoin$group)
pdf("output/beeswarm_plot.pdf", width = 20, height = 10)
plt3 <- ggplot(datjoin, aes(group, yg_n_recount)) + 
  #geom_beeswarm(dodge.width = 0.7, cex = 20) +
  #geom_quasirandom(dodge.width = 0.9, cex = 2, alpha = 0.4, method = "smiley") +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.1)) +
  scale_color_manual(c("gray50", "firebrick")) +
  theme_pubr()
plt3
dev.off()

# Now compare splicease sub-clades

