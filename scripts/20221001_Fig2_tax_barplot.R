# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer", "DECIPHER", "rentrez", "ggpubr", "colorspace",
               "Biostrings", "treeio", "ggtree", "scales", "Biostrings", "readxl", "pals", "svglite")

# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")

# Keep only the putative spliceases
splics <- readLines("data/20210601_spliceases_cut_from_dendro.txt")
seqselect <- tr$tip.label[tr$tip.label %in% splics]
querselect <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
                          paste0(word(seqselect, 1, sep = "_"),
                                 "_", word(seqselect, 2, sep = "_")),
                        TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# acc_phylum_class_order_family_genus_species_gene_name
fixselect <- gsub("WP_|NP_|YP_", "", seqselect)
fixtibb <- tibble(fixselect)
taxdf <- separate(fixtibb, col = fixselect, into = paste0("tax_", 1:8), remove = F, sep = "_", extra = "merge", fill = "right") 
tibbdf <- taxdf %>%
  dplyr::mutate(tax_bar = case_when(grepl("roteobacteria", tax_2) ~ tax_3,
                                    tax_2 %in% c("Unknown") ~ tax_6,
                          TRUE ~ tax_2))

unclass <- tibbdf %>%
  dplyr::filter(tax_6 %in% c("bacterium", "Bacteria", "",
                             "unclassified", "uncultured",
                             "candidate", "Candidatus")) %>%
  pull(tax_1)

quers_to_pull <- querselect[grep(paste0(unclass, collapse = "|"), querselect)]
fetchr <- entrez_fetch(quers_to_pull, db = "protein", rettype = "fasta")
write(fetchr, "~/Documents/github/splicase_analysis/data/unclassified_spliceases.fa")

# Read in and parse seqs
parsr <- readAAStringSet("~/Documents/github/splicase_analysis/data/unclassified_spliceases.fa")
parstibb <- tibble(names(parsr)) 
parstibb

parsdf <- separate(parstibb, col = 1, into = paste0("nam_", 1:8), remove = F, sep = " ", extra = "merge", fill = "right") %>%
  dplyr::mutate(quer_acc = case_when(grepl("WP_|NP_|YP_", nam_1) ~ 
                                       word(nam_1, 2, sep = "_"),
                                     TRUE ~ word(nam_1, 1, sep = "_")))

parsdf$quer_acc <- gsub("\\.1", "", parsdf$quer_acc)
parsdf$quer_acc
tibbdf$tax_1

# Merge data
mergdf <- tibbdf %>%
  left_join(parsdf, by = c("tax_1" = "quer_acc")) %>%
  dplyr::mutate(tax_fill = case_when(tax_6 %in% 
                             c("bacterium", "Bacteria", "",
                             "unclassified", "uncultured",
                             "candidate", "Candidatus") ~ NA_character_,
                             TRUE ~ tax_bar))

# Entrez fetch for the sequences without taxa
tofetch <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
               paste0(word(seqselect, 1, sep = "_"),
                      "_", word(seqselect, 2, sep = "_")),
             TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# Read in the manually-edited taxonomy data
taxdat <- read_excel("data/raw_taxdf_manually_edited.xlsx") %>%
  dplyr::mutate(protein = "Splicease") 

sort(table(taxdat$tax_fill), decreasing = T)
length(table(taxdat$tax_fill))

pal1 <- c(unname(stepped2()), unname(stepped()))
pal2 <- c(tol(), kelly())
pal3 <- darken(pal2, amount = 0.2)
pal.bands(kelly())
pal.bands(tol())
palrem <- c(pal3[3], pal3[4], pal3[7],
          pal3[9], pal3[13], pal3[18], pal3[19], 
          pal3[20], pal3[22])
pal4 <- pal3[!pal3 %in% palrem]

pdf("output/taxa_barplot_horizontal_change_ticks.pdf", width = 7)
ggplot(taxdat, aes(x = protein, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "count") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 22),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = 10) +
  guides(fill=guide_legend(ncol=2)) +
  coord_flip()
dev.off()

sumdat <- taxdat %>%
  group_by(tax_fill) %>%
  add_count(tax_fill) %>%
  dplyr::mutate(perc = n/nrow(.) * 100) %>%
  dplyr::distinct(tax_fill, .keep_all = T)

pdf("output/taxa_barplot_percentage_horizontal_stacked.pdf", width = 10, height = 4)
ggplot(sumdat, aes(x = protein, y = perc, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20), # Left margin
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill=guide_legend(ncol=5)) +
  coord_flip()
dev.off()

colordf <- data.frame(color = pal4,
taxa = sort(unique(sumdat$tax_fill)))
write_csv(colordf, 'data/palette_taxa.csv')
pal4[pal4 == "#A94E35"] <- "gray80"

colnames(sumdat)[colnames(sumdat) == "perc"] <- "Relative abundance (%)"
pdf("output/taxa_barplot_percentage_horizontal_wide_change_ticks.pdf", width = 17, height = 2.5)
ggplot(sumdat, aes(x = protein, y = `Relative abundance (%)`, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.x = element_text( size = 15),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0, 100, by = 10)) +
  guides(fill=guide_legend(ncol=3)) +
  coord_flip()
dev.off()

svglite("output/horizontal_barplot_svg_export.svg", width = 17, height = 2.5)
ggplot(sumdat, aes(x = protein, y = `Relative abundance (%)`, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.x = element_text( size = 15),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0, 100, by = 10)) +
  guides(fill=guide_legend(ncol=3)) +
  coord_flip()
dev.off()

