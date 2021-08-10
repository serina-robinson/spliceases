# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer", "DECIPHER", "rentrez", "ggpubr", "colorspace",
               "Biostrings", "treeio", "ggtree", "scales", "Biostrings", "readxl", "pals")

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
table(tibbdf$tax_bar)

write_csv(tibbdf, "data/raw_taxdf.csv")

unclass <- tibbdf %>%
  dplyr::filter(tax_6 %in% c("bacterium", "Bacteria", "",
                             "unclassified", "uncultured",
                             "candidate", "Candidatus")) %>%
  pull(tax_1)

quers_to_pull <- querselect[grep(paste0(unclass, collapse = "|"), querselect)]
quers_to_pull
#fetchr <- entrez_fetch(quers_to_pull, db = "protein", rettype = "fasta")
#write(fetchr, "~/Documents/github/splicase_analysis/data/unclassified_spliceases.fa")

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

# Merge with the previous df 
parsdf$nam_1
mergdf <- tibbdf %>%
  left_join(parsdf, by = c("tax_1" = "quer_acc")) %>%
  dplyr::mutate(tax_fill = case_when(tax_6 %in% 
                             c("bacterium", "Bacteria", "",
                             "unclassified", "uncultured",
                             "candidate", "Candidatus") ~ NA_character_,
                             TRUE ~ tax_bar))

write_csv(mergdf, "data/raw_taxdf_to_name.csv")


# Entrez fetch for the sequences without taxa
tofetch <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
               paste0(word(seqselect, 1, sep = "_"),
                      "_", word(seqselect, 2, sep = "_")),
             TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# Filter data for spliceases
# gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv") %>%
#   dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) %>%
#   dplyr::filter(row_id %in% c(3, 4, 5, 7, 8, 9))

# Read in the manually-edited taxonomy data
taxdat <- read_excel("data/raw_taxdf_manually_edited.xlsx") %>%
  dplyr::mutate(protein = "Splicease") #%>%
  #dplyr::filter(Rel_abundance = tax_fill)

sort(table(taxdat$tax_fill), decreasing = T)
length(table(taxdat$tax_fill))

pal1 <- c(unname(stepped2()), unname(stepped()))
pal2 <- c(tol(), kelly())
pal3 <- darken(pal2, amount = 0.2)
pal.bands(kelly())
pal.bands(tol())
palrem <- c(pal3[3], pal3[4], pal3[7],
          pal3[9], pal3[13], pal3[18], pal3[19], 
          pal3[20], pal3[22]) #"#F2F3F4", "#F3C300",
        #  "#875692", "#6E8F01", "#A1CAF1", "#BE0032")
pal4 <- pal3[!pal3 %in% palrem]


pdf("output/taxa_barplot_2col.pdf", width = 7)
ggplot(taxdat, aes(x = protein, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "count") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill=guide_legend(ncol=2))
  #scale_y_continuous(limits = c(0,1000), breaks = c(seq(from = 0, to = 1000, by = 100)))
  # theme(axis.text.x = element_text(angle = 75, hjust = 1),
  #       #   axis.text.y = element_blank(),
  #       plot.margin = margin(2, 2, 2, 2, "cm"))
dev.off()


sumdat <- taxdat %>%
  group_by(tax_fill) %>%
  add_count(tax_fill) %>%
  dplyr::mutate(perc = n/nrow(.) * 100) %>%
  dplyr::distinct(tax_fill, .keep_all = T)
sum(sumdat$perc)

pdf("output/taxa_barplot_percentage.pdf", width = 6.5)
ggplot(sumdat, aes(x = protein, y = perc, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = pal4) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill=guide_legend(ncol=2))
#scale_y_continuous(limits = c(0,1000), breaks = c(seq(from = 0, to = 1000, by = 100)))
# theme(axis.text.x = element_text(angle = 75, hjust = 1),
#       #   axis.text.y = element_blank(),
#       plot.margin = margin(2, 2, 2, 2, "cm"))
dev.off()
  