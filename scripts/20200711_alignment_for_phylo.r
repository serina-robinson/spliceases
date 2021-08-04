# Install packages
pacman::p_load("data.table", "tidyverse", "Biostrings", "treeio")

# New BLAST results
sqs <- readAAStringSet("data/blast_results/1625_blast_hits_seqs_pulled.fasta.txt")

# Find the taxonomic assignments of the sequences
parsed_tax <- fread("data/taxonomy/parsed_tax.tsv", data.table = F) %>%
  janitor::clean_names() %>%
  dplyr::distinct(genus, .keep_all=TRUE)

tax_df <- data.frame(nam = names(sqs),
                     aa_seq = sqs) %>%
  dplyr::mutate(acc = stringr::word(nam, 1, sep = "\\.1")) %>%
  dplyr::mutate(specnam = gsub("\\]", " ", stringr::word(nam, 2, sep = "\\["))) %>%
  dplyr::mutate(genus = stringr::word(specnam, 1, sep = " "))
tax_df$genus

# Merge with the taxonomy
ptax <- tax_df %>%
  dplyr::left_join(., parsed_tax, by = "genus")

no_tax <- ptax$genus[is.na(ptax$phylum)]
no_tax

# Define a naming system which is as follows
# acc_phylum_class_order_family_genus_species_gene_name

# Fix names
tree_tax <- ptax %>%
  dplyr::mutate(sqnams = paste(acc, phylum, class, order, family, genus, "BLAST", sep = "_"))

tree_tax$sqnams <- gsub("-", "_", tree_tax$sqnams)
tree_tax$sqnams <- gsub("\\.", "_", tree_tax$sqnams)
tree_tax$sqnams <- gsub(" ", "_", tree_tax$sqnams)
tree_tax$sqnams <- gsub("_NA_", "_Unknown_", tree_tax$sqnams)
tree_tax$sqnams <- gsub("_No_class_", "_Unknown_", tree_tax$sqnams)
tree_tax$sqnams

# To write to file
new <- AAStringSet(tree_tax$aa_seq)
names(new) <- tree_tax$sqnams
length(new)
writeXStringSet(new, "output/1625_BLAST_hits_seqs_unaligned.fa")

# Read in the old unaligned sequence set
old <- readAAStringSet("data/seqs_for_phylogeny/1796_combined_sqs_for_phylo_v1_fixnams.fa")
comb <- AAStringSet(c(new, old))
dedup <- comb[!duplicated(comb)]
summary(width(dedup)) # min 157, max 782
writeXStringSet(dedup, "data/seqs_for_phylogeny/3418_expanded_blast_seqs_for_phylo_names_fixed.fa")

# Align using MUSCLE
# muscle -in 3418_expanded_blast_seqs_for_phylo_names_fixed.fa -out 3418_expanded_blast_seqs_for_phylo_names_fixed_2iters.afa -maxiters 2

# Trim alignment using trimal
# trimal -in 3418_expanded_blast_seqs_for_phylo_names_fixed_2iters.afa -gt 0.1 -htmlout 3418_aligned.html > 3418_expanded_blast_seqs_for_phylo_v2_aligned_2iters_trimmed.fasta

# Approx. maximum-likelihood with FastTree
# fasttree <3418_expanded_blast_seqs_for_phylo_v2_aligned_2iters_trimmed.fasta> 3418_fasttree_gt_10.nwk

aln_raw <- readAAStringSet("data/seqs_for_phylogeny/3418_expanded_blast_seqs_for_phylo_names_fixed_2iters.afa")
aln <- readAAStringSet("data/seqs_for_phylogeny/3418_expanded_blast_seqs_for_phylo_v2_aligned_2iters_trimmed.fasta")

# Final phylogenetic tree
tr <- read.tree("data/seqs_for_phylogeny/3418_fasttree_gt_10.nwk")
