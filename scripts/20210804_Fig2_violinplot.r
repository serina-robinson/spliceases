# Install packages
pacman::p_load("tidyverse", "ggstatsplot", "ggsignif",
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
  left_join(., grps, by = "query_acc") %>%
  dplyr::mutate(is_splicease = ifelse(query_acc %in% query_rem, "Spliceases", "Other rSAM-SPASM")) %>%
  dplyr::mutate(aa_len = nchar(aa)) %>%
  dplyr::mutate(uniq_id = paste0(protein_acc, "_", variable)) %>%
  dplyr::filter(!duplicated(uniq_id)) %>%
  dplyr::filter(row_id %in% c(3,4,5,7,8,9)) %>%
  dplyr::group_by(query) %>%
  dplyr::mutate(yg_density = (sum(yg_n_recount)/sum(aa_len))) %>%
  ungroup()

# Clean up for plotting
datjoin$group[is.na(datjoin$group)] <- 11
datjoin$group <- as.factor(datjoin$group)

pdf("output/splicease_non_splicease_plot_density_6genes.pdf", width = 4, height = 3.5)
plt3 <- ggplot(datjoin, aes(is_splicease, yg_density)) +
  geom_violin(alpha = 0.4, aes(color = is_splicease, fill = is_splicease)) +
  geom_point(position = position_jitter(seed = 1, width = 0.1), alpha = 0.5,
             aes(color = is_splicease), size = 1) +
  scale_color_manual(values = c("gray50", "#E41A1C")) +
  scale_fill_manual(values = c("gray50", "#E41A1C")) +
  theme_pubr() +
  ylab("XYG frequency") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  geom_signif(
    comparisons = list(c("Spliceases", "Other rSAM-SPASM")),
    color = "black",
  map_signif_level = function(p) sprintf("*** p = %.2g", p))
plt3
dev.off()


# View the distribution of the dataset for statistical testing
hist(log(datjoin$yg_n_recount[datjoin$is_splicease == "Splicease"]))
hist(log(datjoin$yg_density[datjoin$is_splicease == "Splicease"]))
hist(log(datjoin$yg_n_recount[datjoin$is_splicease == "Other rSAM-SPASM"]))
hist(log(datjoin$yg_density[datjoin$is_splicease == "Other rSAM-SPASM"]))

# Welch's t-test
dens_diff <- t.test(log(datjoin$yg_density[datjoin$is_splicease == "Splicease"]),
                    log(datjoin$yg_density[datjoin$is_splicease == "Other rSAM-SPASM"]),
                    alternative = "two.sided",
                    var.equal = FALSE)
dens_diff


# unpaired two-samples Wilcoxon test 
count_diff <- wilcox.test(datjoin$yg_n_recount[datjoin$is_splicease == "Other rSAM-SPASM"],
                          datjoin$yg_n_recount[datjoin$is_splicease == "Splicease"],
                          alternative = "two.sided",
                          exact = FALSE)
count_diff 

# unpaired two-samples Wilcoxon test 
res <- wilcox.test(splicdat$yg_density, 
                   otherdat$yg_density,
                   alternative = "two.sided", # also greater 
                   exact = FALSE)



