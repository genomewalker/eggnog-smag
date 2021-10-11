library(tidyverse)
library(RSQLite)
library(maditr)
library(pbmcapply)
library(tidytext)
library(tidygramr)
library(stringdist)

unixtools::set.tempdir(path.expand("/vol/cloud/antonio/tmp"))
Sys.setenv(RENV_PATHS_CACHE = "/vol/cloud/antonio/tmp/renv/cache")

db <- "/dev/shm/eggnog-smags.sqlite"
con <- RSQLite::dbConnect(RSQLite::SQLite(), db)


eggdb <- "data/eggnog.db"
con1 <- RSQLite::dbConnect(RSQLite::SQLite(), eggdb)

results_dir <- "data/01_Results_EggNog_annotationV1"
emapper_files_pattern <- ".metaeuk.eggnog.annotation.emapper.annotations"
emapper_files_pattern_filt <- paste0(emapper_files_pattern, "-filt.tsv")
emapper_files <- list.files(path = results_dir, pattern = emapper_files_pattern, full.names = TRUE)

taxmap <- read_tsv("data/eggnog-db/emapper-v2.taxmap.tsv",
                   col_names = c("taxid", "tax_name"),
                   col_types = list("c", "c"), comment = "#")

funcat <- fread("data/eggnog-db/eggnog-funcat.tsv",) %>%
  setNames(c("COG_Functional_Category", "cog_category", "cog_description"))



read_emmaper_files <- function(file, pattern = pattern){
  col_names <- c("query_name","seed_eggNOG_ortholog","seed_ortholog_evalue",
                 "seed_ortholog_score","Predicted_taxonomic_group","Predicted_protein_name",
                 "Gene_Ontology_terms","EC_number","KEGG_ko","KEGG_Pathway","KEGG_Module",
                 "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","tax_scope",
                 "eggNOG_OGs","bestOG","COG_Functional_Category","eggNOG_free_text_description")
  fread(file, col.names = col_names, quote = "", sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(pattern = pattern, replacement = "", x = basename(file))) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(eggNOG_OGs, into = c("og", "taxid"), sep = "@", remove = FALSE) %>%
    left_join(taxmap, by = "taxid") %>%
    filter(Predicted_taxonomic_group == tax_name) %>%
    distinct() %>%
    left_join(annotations) %>%
    mutate(bestOG = og)
  # if (dbExistsTable(con, "emapper_results") == TRUE){
  #   DBI::dbWriteTable(conn = con,
  #                     name = "emapper_results",
  #                     value = res,
  #                     append= TRUE)
  # }else{
  #   DBI::dbWriteTable(conn = con,
  #                     name = "emapper_results",
  #                     value = res)
  # }
}

read_emmaper_files_lower <- function(file, pattern = pattern){
  col_names <- c("query_name","seed_eggNOG_ortholog","seed_ortholog_evalue",
                 "seed_ortholog_score","Predicted_taxonomic_group","Predicted_protein_name",
                 "Gene_Ontology_terms","EC_number","KEGG_ko","KEGG_Pathway","KEGG_Module",
                 "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","tax_scope",
                 "eggNOG_OGs","bestOG","COG_Functional_Category","eggNOG_free_text_description")
  fread(file, col.names = col_names, quote = "", sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(pattern = pattern, replacement = "", x = basename(file))) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(eggNOG_OGs, into = c("og", "taxid"), sep = "@", remove = FALSE) %>%
    left_join(taxmap, by = "taxid") %>%
    group_by(sMAG, query_name) %>%
    #slice(which.min(taxid)) %>% View()
    filter(taxid == min(taxid)) %>%
    ungroup() %>%
    #filter(Predicted_taxonomic_group == tax_name) %>%
    distinct() %>%
    select(sMAG, query_name, og, taxid) %>%
    distinct() %>%
    mutate(bestOG = og) %>%
    inner_join(annotations %>% rename(taxid = level)) %>%
    mutate(description = ifelse(description == "", "Not available", description))
  # if (dbExistsTable(con, "emapper_results") == TRUE){
  #   DBI::dbWriteTable(conn = con,
  #                     name = "emapper_results",
  #                     value = res,
  #                     append= TRUE)
  # }else{
  #   DBI::dbWriteTable(conn = con,
  #                     name = "emapper_results",
  #                     value = res)
  # }
}


annotations <- tbl(con1, "og") %>% collect()

if (dbExistsTable(con, "emapper_results") == TRUE){
  dbRemoveTable(con, "emapper_results")
}


emapper_results <- pbmclapply(emapper_files, read_emmaper_files, pattern = emapper_files_pattern, mc.cores = 28)
emapper_results_lower <- pbmclapply(emapper_files, read_emmaper_files_lower, pattern = emapper_files_pattern, mc.cores = 16)

emapper_results <- emapper_results %>% bind_rows()
emapper_results_lower <- emapper_results_lower %>% bind_rows()

smags <- emapper_results$sMAG %>% unique()
smags_lower <- emapper_results_lower$sMAG %>% unique()


pbmclapply(smags_lower, function(X){
  fname <- file.path("results/lowestOG/01_Results_EggNog_annotationV1-filt", paste0(X, emapper_files_pattern_filt))
  emapper_results_lower %>%
    filter(sMAG == X) %>%
    write_tsv(path = fname)
}, mc.cores = 8)

DBI::dbWriteTable(conn = con,
                  name = "emapper_results",
                  value = emapper_results,
                  append= TRUE)


emapper_results_agg <- emapper_results_lower %>%
  inner_join(text) %>%
  rename(description_clean = lines) %>%
  select(sMAG, description_clean) %>%
  group_by(sMAG, description_clean) %>%
  count(sort = T) %>%
  ungroup()

DBI::dbWriteTable(conn = con,
                  name = "emapper_results_agg",
                  value = emapper_results_agg,
                  append= TRUE)

write_tsv(emapper_results_agg, "results/lowestOG/01_Results_EggNog_annotationV1-filt_agg-clean.tsv.gz")

emapper_results_descs <- emapper_results %>%
  select(eggNOG_free_text_description, description) %>%
  distinct()

emapper_results_descs %>%
  filter(eggNOG_free_text_description != description) %>%
  nrow()

text <- emapper_results_lower %>%
  select(description) %>%
  distinct() %>%
  mutate(text = gsub("_", " ", description)) %>%
  unnest_tokens(lines, text, token = "lines", drop = TRUE)

sep <-  paste(rep("#", 5), collapse = "")
text_comb <- text %>%
  group_by(lines) %>%
  summarise(comb_descriptions = paste(description, collapse = sep))

derep <- text_comb %>%
  filter(grepl(sep, comb_descriptions))
rr <- text_comb %>%
  filter(!(grepl(sep, comb_descriptions)))

rr_suffix <- rr %>%
  filter(grepl("[- ][A-Za-z]$", lines)) %>%
  bind_rows(rr %>%
              filter(grepl("[- ][0-9]\\d*$", lines))) %>%
  bind_rows(rr %>%
              filter(grepl(" i*$", lines))) %>%
  bind_rows(rr %>%
              filter(grepl("duf", lines)))

rr_filt <- rr %>%
  filter(!lines %in% rr_suffix$lines) %>% arrange(lines)


derep %>%
  select(comb_descriptions) %>%
  separate_rows(comb_descriptions, sep = sep)
#
# words <- text_comb %>%
#   unnest_tokens(word, lines, token = "regex",
#                 to_lower = TRUE,
#                 pattern = "\\s+",
#                 drop = FALSE)
#
# one_grams <- words %>%
#   group_by(comb_descriptions) %>%
#   count() %>%
#   filter(n < 2)
#
# text_1 <- text %>%
#   filter(!(lines %in% c(one_grams$comb_descriptions, derep$lines)))
#
# bigrams <- text_1 %>%
#   unnest_tokens(bigram, lines, token = "ngrams", n = 2, drop = FALSE)
#
# bigrams_unique <- bigrams %>%
#   group_by(lines) %>%
#   count() %>%
#   filter(n < 2)
#
#
# text_2 <- text_1 %>%
#   filter(!(lines %in% c(bigrams_unique$lines)))
#
# trigrams <- text_2 %>%
#   unnest_tokens(trigram, lines, token = "ngrams", n = 3, drop = FALSE)
#
# trigrams_unique <- trigrams %>%
#   group_by(lines) %>%
#   count() %>%
#   filter(n < 2)

d <- rr_filt %>%
  ungroup() %>%
  mutate(text = lines) %>%
  mutate(text = gsub("\\s+", "", lines, perl = TRUE)) %>%
  clean_corpus(rm_special = TRUE)

da_dist <- stringdist::stringdistmatrix(d$text %>% unique(), useNames = TRUE, method = "jaccard", q = 7, nthread = 28)

da_dist_m <- as.matrix(da_dist)
da_dist_m[upper.tri(da_dist_m)] <- NA
da_dist_m[da_dist_m > 0.2] <- NA

idx <- which(!is.na(da_dist_m), arr.ind = T)
df <- as.data.frame(da_dist_m)
cols <- names(df)
#rows <- rownames(df)

#cols %>% enframe() %>% filter(value == "X1phosphatidylinositolbinding")

idx <- as.data.table(idx) %>% take_if(row != col)

get_idx <- function(X){
  tibble(row = X[[1]], col = X[[2]], item1=cols[X[[1]]], item2 = cols[X[[2]]], dist=df[X[[1]],X[[2]]])
}

pb <- txtProgressBar(min = 0, max = nrow(idx), style = 3)
idx <- idx[, {
  setTxtProgressBar(pb, .GRP);
  get_idx(.SD);
}, by =  1:nrow(idx), .SDcols = c("row", "col")]
close(pb)

idx %>%
  as_tibble() %>%
  arrange(dist) %>% View()


test1 <- test %>% as_tibble %>%
  filter(item1 != variable, !is.na(value)) %>% filter(value < 0.1) %>% arrange(desc(value))

test1 %>%
  inner_join(d %>% rename(item1 = text) %>% select(-n)) %>%
  rename(linesA = lines) %>%
  inner_join(d %>% rename(variable = text) %>% select(-n))

d %>% filter(text %in% (c(test1$item1, test1$variable) %>% unique()))

d %>% filter(text == "V1269")
