# Install packages
pacman::p_load("data.table", "tidyverse", "RColorBrewer", "DECIPHER",
               "Biostrings", "treeio", "ggtree", "scales", "Biostrings")

# Read in the sequences
testseqs <- readAAStringSet('data/seqs_for_stats/All_flanking_genes_5_7_for_Tom.fasta')
names(testseqs)[grep("archa", names(testseqs))] # archaea are included
length(testseqs) # 794 sequences

# Read in the consensus tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")

# Keep only the putative spliceases
torem <- readLines("data/20210601_spliceases_cut_from_dendro.txt")
torem
tokeep <- tr$tip.label[!tr$tip.label %in% torem]

set.seed(1234)
#rand <- sample(x = length(tokeep), size = 794, replace = F)
#seqselect <- tokeep[rand]
seqselect <- tokeep
tokeep
querselect <- case_when(grepl("WP_|NP_|YP_", seqselect) ~ 
                        paste0(word(seqselect, 1, sep = "_"),
                        "_", word(seqselect, 2, sep = "_")),
                        TRUE ~ paste0(word(seqselect, 1, sep = "_")))


gbardat <- read_csv("data/XYG_count_data/combined_YG_count_data.csv") %>%
  dplyr::filter(grepl(paste0(querselect, collapse = "|"), query_acc)) %>%
  #dplyr::filter(query_acc %in% querselect) %>%
  dplyr::filter(row_id %in% c(5,7)) %>%
  distinct(aa, .keep_all = T)

ctrl_aa <- AAStringSet(gbardat$aa)
names(ctrl_aa) <- paste0(gbardat$row_id, "_", gbardat$query_acc, "_", gbardat$protein_acc, 
                                               "_", gbardat$genus_species, "_", gbardat$name1)
# writeXStringSet(ctrl_aa, "data/533_control_rSAM_SPASM_seqs.fasta")

# MEME suite
# https://meme-suite.org/meme/info/status?service=AME&id=appAME_5.3.31628507280967550018493
# ame --verbose 1 --oc . --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 All_flanking_genes_5_7_for_Tom.fasta motifs.meme

# rank	motif_DB	motif_ID	motif_alt_ID	consensus	p-value	adj_p-value	E-value	tests	FASTA_max	pos	neg	PWM_min	TP	%TP	FP	%FP
# 1	motifs.meme	1	PYG	PYG	1.03e-11	8.15e-9	8.15e-9	793	379	379	415	1.03	77	20.32	20	4.82

# AME (Analysis of Motif Enrichment): Version 5.3.3 compiled on Feb 21 2021 at 14:51:06
# The format of this file is described at https://meme-suite.org/meme/doc/ame-output-format.html.
# ame --verbose 1 --oc . --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 All_flanking_genes_5_7_for_Tom.fasta motifs.meme



