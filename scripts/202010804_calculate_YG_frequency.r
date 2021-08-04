# Install packages
pacman::p_load("ggplot2", "tidyverse", "Biostrings", "ggpubr",
               "rentrez", "mltools", "data.table", "caret")

# Read in the characterized seq dataset
char1 <- read_csv("data/rodeo2_output/characterized_main_co_occur.csv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(acc_cln = word(protein_acc, sep = "\\.1", 1)) %>%
  dplyr::group_by(query) %>%
  add_tally(., name = "nb_size") %>%
  dplyr::mutate(row_id = row_number()) %>%
  ungroup() %>%
  select(-contains("x"))

# Read in the characterized sequences
sqs1 <- readAAStringSet("data/characterized_homologs/characterized_seqs_nbs.fasta")

# Match up the two datasets
sqdf1 <- data.frame(nams = names(sqs1),
                   aa = sqs1) %>%
  dplyr::mutate(acc_cln = word(nams, 1, sep = "\\.1")) %>%
  distinct(acc_cln, .keep_all = T)

joindf1 <- char1 %>%
  left_join(., sqdf1, by = "acc_cln") 

# Read in the BLAST dataset
char2 <- read_csv("data/rodeo2_output/1703_BLAST_main_co_occur.csv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(acc_cln = word(protein_acc, sep = "\\.1", 1)) %>%
  dplyr::group_by(query) %>%
  add_tally(., name = "nb_size") %>%
  dplyr::mutate(row_id = row_number()) %>%
  ungroup()

# Read in the BLAST sequences
sqfils <- list.files("data/blast_results/blast_gene_neighborhoods", full.names = T)

aa_comb <- AAStringSet()
for(i in 1:length(sqfils)){
  aa_comb <- AAStringSet(c(aa_comb, readAAStringSet(sqfils[i])))
}

sqs2 <- aa_comb[!duplicated(aa_comb)]

# Match up the two datasets
sqdf2 <- data.frame(nams = names(sqs2),
                   aa = sqs2) %>%
  dplyr::mutate(acc_cln = word(nams, sep = "\\.1", 1)) %>%
  distinct(acc_cln, .keep_all = T)

joindf2 <- char2 %>%
  left_join(., sqdf2, by = "acc_cln") 

# Read in the key dataset
char3 <- read_csv("data/rodeo2_output/58_SPASM_key_main_co_occur.csv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(acc_cln = word(protein_acc, sep = "\\.1", 1)) %>%
  dplyr::group_by(query) %>%
  add_tally(., name = "nb_size") %>%
  dplyr::mutate(row_id = row_number()) %>%
  ungroup()

# Read in the key sequences
sqfils3 <- list.files("data/spasm_key/key_gene_neighborhoods/", full.names = T)
sqfils3
aa_comb <- AAStringSet()
for(i in 1:length(sqfils3)){
  aa_comb <- AAStringSet(c(aa_comb, readAAStringSet(sqfils3[i])))
}

sqs3 <- aa_comb[!duplicated(aa_comb)]

# Match up the two datasets
sqdf3 <- data.frame(nams = names(sqs3),
                   aa = sqs3) %>%
  dplyr::mutate(acc_cln = word(nams, 1, sep = "\\.1")) %>%
  distinct(acc_cln, .keep_all = T)

joindf3 <- char3 %>%
  left_join(., sqdf3, by = "acc_cln") 


# Read in the PFAM dataset
char4 <- read_csv("data/rodeo2_output/1886_main_co_occur.csv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(acc_cln = word(protein_acc, sep = "\\.1", 1)) %>%
  dplyr::group_by(query) %>%
  add_tally(., name = "nb_size") %>%
  dplyr::mutate(row_id = row_number()) %>%
  ungroup()

# Read in the sequences
sqfils <- list.files("data/1886_refseq_gene_neighborhoods", full.names = T)

aa_comb <- AAStringSet()
for(i in 1:length(sqfils)){
  aa_comb <- AAStringSet(c(aa_comb, readAAStringSet(sqfils[i])))
}

sqs4 <- aa_comb[!duplicated(aa_comb)]

# Match up the two datasets
sqdf4 <- data.frame(nams = names(sqs4),
                   aa = sqs4) %>%
  dplyr::mutate(acc_cln = word(nams, 1, sep = "\\.1")) %>%
  distinct(acc_cln, .keep_all = T)

joindf4 <- char4 %>%
  left_join(., sqdf4, by = "acc_cln") 

# Collect YA|YG co-occurence
ygadf <- joindf1 %>%
  bind_rows(joindf2, joindf3, joindf4) %>%
  dplyr::mutate(yg = gregexpr(pattern = ".YG", aa)) %>%
  dplyr::mutate(ya = gregexpr(pattern = ".YA", aa)) %>%
  separate(yg, into = paste0("yg_", 1:20), remove = F, sep = ",", extra = "merge", fill = "right") %>%
  separate(ya, into = paste0("ya_", 1:20), remove = F, sep = ",", extra = "merge", fill = "right") %>%
  mutate_at(vars(c(paste0("yg_", 1:20), paste0("ya_", 1:20))), ~ str_replace(., "c\\(|\\)", "")) %>%
  mutate_at(vars(c(paste0("yg_", 1:20), paste0("ya_", 1:20))), ~ str_replace(., "-1", NA_character_)) %>%
  dplyr::mutate(query_acc = stringr::word(query, sep = "\\.1", 1)) %>%
  mutate_at(vars(c(paste0("yg_", 1:20), paste0("ya_", 1:20))), ~ as.numeric(.)) %>%
  mutate_at(vars(c(paste0("yg_", 1:20), paste0("ya_", 1:20))), ~ substr(x = aa, start = ., stop = .)) %>%
  dplyr::select(-ya, -yg)

# Now need to summarize YGs somehow...so there is only one row per query...
# First look at 'YG' motif occurence in SAMs
ygasam <- ygadf %>%
  dplyr::filter(protein_acc == query) # Radical SAM YGs

# Remove radical SAMs
yga_total <- ygadf %>%
  dplyr::filter(protein_acc != query)

# Convert to long format
ygalong <- yga_total %>%
  dplyr::select(row_id, nb_size, query, genus_species, nucleotide_acc, protein_acc,
                name1, pfam_id1, description1, aa, paste0("yg_", 1:20)) %>% #paste0("ya_", 1:20)) %>% # if want YAs too
  reshape2::melt(., id.vars = c("row_id", "nb_size", "query", "genus_species", "nucleotide_acc", "protein_acc",
                                "name1", "pfam_id1", "description1", "aa")) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(query_acc = stringr::word(query, sep = "\\.1", 1)) %>%
  dplyr::filter(!duplicated(.))

ygabar <- ygalong %>%
  dplyr::distinct(., query, variable, .keep_all =T) %>%
  group_by(query_acc) %>%
  dplyr::add_count() %>%
  dplyr::mutate(xlab = paste0(query_acc, " ", genus_species)) %>%
  dplyr::mutate(n = as.numeric(n)) %>%
  dplyr::mutate(yg_n_recount = str_count(aa, pattern = ".YG")) 

write_csv(ygabar, "data/XYG_count_data/combined_YG_count_data.csv")



