# Install packages
pacman::p_load("ggtree", "stringr", "ggrepel", "dplyr", "scales", "ape", "DECIPHER", "readxl", "viridis",
               "ggtreeExtra", "ggnewscale", "readr", "RColorBrewer", "phytools", "rentrez", "ggpubr")

# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
torem <- readLines("data/20210601_spliceases_cut_from_dendro.txt") 
newtr <- keep.tip(tr, grep(paste0(torem, collapse = "|"), tr$tip.label))
midtr <- midpoint.root(newtr)

pdf("output/rectangular_geom_tiplabs.pdf", width = 50, height = 120)
gdf <- ggtree(midtr, layout = "rectangular", color = "gray60") +
    geom_tiplab(aes(size = 2)) +
    geom_nodelab(aes(label = node, size = 2 , color = "red")) +
    xlim(NA, 10)
gdf 
dev.off()

# Extract the GDL clade
gl <- groupClade(midtr, .node = 1288)
gtfake <- ggtree(gl)
tips <- gtfake$data$label[gtfake$data$isTip & gtfake$data$group == 1]
tips
genera <- word(tips, sep = "_", -2)
table(genera) 


accs <- case_when(grepl("WP_|YP_|NP_", tips) ~ paste0(word(tips, sep = "_", 1), "_", word(tips, sep = "_", 2)),
                  TRUE ~ paste0(word(tips, sep = "_", 1)))
writeLines(accs, "data/20222001_Caulobacter_GDL_accs.txt")

# Now read in the GDL phylogeny only
# Keep only the putative spliceases
splics <- readLines("data/20222001_Caulobacter_GDL_accs.txt")
# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")

# Keep only the putative spliceases
seqselect <- tr$tip.label[grep(paste0(splics, collapse = "|"), tr$tip.label)]

querselect <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
                          paste0(word(seqselect, 1, sep = "_"),
                                 "_", word(seqselect, 2, sep = "_")),
                        TRUE ~ paste0(word(seqselect, 1, sep = "_")))

# acc_phylum_class_order_family_genus_species_gene_name
fixselect <- gsub("WP_|NP_|YP_", "", seqselect)
fixtibb <- tibble(fixselect)
fixtibb
taxdf <- tidyr::separate(fixtibb, col = fixselect, into = paste0("tax_", 1:8), remove = F, sep = "_", extra = "merge", fill = "right") 
tibbdf <- taxdf %>%
  dplyr::mutate(tax_bar = case_when(grepl("roteobacteria", tax_2) ~ tax_3,
                                    tax_2 %in% c("Unknown") ~ tax_6,
                                    TRUE ~ tax_2))
unclass <- tibbdf %>%
  dplyr::filter(tax_6 %in% c("bacterium", "Bacteria", "",
                             "unclassified", "uncultured",
                             "candidate", "Candidatus")) %>%
  pull(tax_1)
unclass
quers_to_pull <- querselect[grep(paste0(unclass, collapse = "|"), querselect)]
quers_to_pull
fetchr <- entrez_fetch(quers_to_pull, db = "protein", rettype = "fasta")
write(fetchr, "~/Documents/github/splicase_analysis/data/unclassified_GDL_substrates.fa")

# Read in and parse seqs
parsr <- readAAStringSet("~/Documents/github/splicase_analysis/data/unclassified_GDL_substrates.fa")
parstibb <- tibble(names(parsr)) 
parstibb

parsdf <- tidyr::separate(parstibb, col = 1, into = paste0("nam_", 1:8), remove = F, sep = " ", extra = "merge", fill = "right") %>%
  dplyr::mutate(quer_acc = case_when(grepl("WP_|NP_|YP_", nam_1) ~ 
                                       word(nam_1, 2, sep = "_"),
                                     TRUE ~ word(nam_1, 1, sep = "_")))
parsdf
parsdf$quer_acc <- gsub("\\.1", "", parsdf$quer_acc)


# Merge with the previous df 
mergdf <- tibbdf %>%
  left_join(parsdf, by = c("tax_1" = "quer_acc")) %>%
  dplyr::mutate(tax_fill = case_when(tax_6 %in% 
                                       c("bacterium", "Bacteria", "",
                                         "unclassified", "uncultured",
                                         "candidate", "Candidatus") ~ NA_character_,
                                     TRUE ~ tax_bar))
write_csv(mergdf, "data/GDL_Caulobacter_raw_taxdf_to_name.csv")

# Read in the manually-edited taxonomy data
taxdat <- read_excel("data/GDL_Caulobacter_raw_taxdf_to_name.xlsx") %>%
  dplyr::mutate(protein = "") %>%
  dplyr::filter(tax_6 != "Acidobacteria")

sort(table(taxdat$tax_fill), decreasing = T)
table(taxdat$tax_fill)

# Read in the palette 
pal_csv <- read_csv("data/palette_taxa.csv")
pal_csv[1:20,]
pal_trim <- pal_csv %>%
  dplyr::filter(taxa %in% unique(taxdat$tax_fill)) %>%
  pull(color)

pdf("output/GDL_taxa_barplot_horizontal.pdf", width = 11, height = 8)
ggplot(taxdat, aes(x = protein, fill = tax_fill)) +
  geom_bar(position = "stack", stat = "count") +
  scale_fill_manual(values = pal_trim) +
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
  guides(fill=guide_legend(ncol=2)) +
  coord_flip()
dev.off()

table(taxdat$tax_fill)


gdl <- data.frame(table(taxdat$tax_fill))
colnames(gdl) <- c("Taxa", "Count")
pdf("output/GDL_donut.pdf")
donut <- ggdonutchart(data = gdl,
                      x = "Count",
                      fill = "Taxa",
                      lab.font = c(8, "plain", "black"),
                      color = "white",
                      palette = pal_trim)
donut
dev.off()



gdl <- data.frame(table(taxdat$tax_fill))

nrow(gdl)
colnames(gdl) <- c("Taxa", "Count")
pdf("output/GDL_donut_phylum.pdf", width = 10)
donut <- ggdonutchart(data = gdl,
                      x = "Count",
                      label = rep("", nrow(gdl)),
                      fill = "Taxa",
                      lab.font = c(8, "plain", "black"),
                      color = "white",
                      palette = pal_trim)
donut
dev.off()



gdl <- data.frame(table(taxdat$tax_6))
pal2 <- inferno(41)

colnames(gdl) <- c("Taxa", "Count")
pdf("GDL_donut_genus.pdf", width = 9)
donut <- ggdonutchart(data = gdl,
                      x = "Count",
                     label = rep("", nrow(gdl)),
                      fill = "Taxa",
                      lab.font = c(8, "plain", "white"),
                      color = "white",
                      palette = pal2)
donut
dev.off()

labs <- paste0(round((gdl$Count/sum(gdl$Count) * 100), 2), "%")
labs
pdf("output/GDL_donut_genus_unlabeled.pdf", width = 9)
donut <- ggdonutchart(data = gdl,
                      x = "Count",
                      label = rep("", nrow(gdl)),
                      fill = "Taxa",
                      lab.font = c(8, "plain", "white"),
                      color = "white",
                      palette = pal2)
donut
dev.off()

pdf("output/GDL_genus_barplot_horizontal.pdf", width = 18, height = 2.5)
ggplot(taxdat, aes(x = protein, fill = tax_6)) +
  geom_bar(position = "stack", stat = "count") +
  scale_fill_manual(values = pal2) +
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
  guides(fill=guide_legend(ncol=5)) +
  coord_flip()
dev.off()



