library(tidyverse)
library(data.table)
# Read and process the emmaper results - RAW
read_emmaper_files <- function(file, pattern = pattern) {
  col_names <-
    c(
      "query_name",
      "seed_eggNOG_ortholog",
      "seed_ortholog_evalue",
      "seed_ortholog_score",
      "Predicted_taxonomic_group",
      "Predicted_protein_name",
      "Gene_Ontology_terms",
      "EC_number",
      "KEGG_ko",
      "KEGG_Pathway",
      "KEGG_Module",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "tax_scope",
      "eggNOG_OGs",
      "bestOG",
      "COG_Functional_Category",
      "eggNOG_free_text_description"
    )
  fread(file,
        col.names = col_names,
        quote = "",
        sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(
      pattern = pattern,
      replacement = "",
      x = basename(file)
    )) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(
      eggNOG_OGs,
      into = c("og", "taxid"),
      sep = "@",
      remove = FALSE
    ) %>%
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

# Read and process the emmaper results - LOWEST OG
read_emmaper_files_lower <- function(file, pattern = pattern) {
  col_names <-
    c(
      "query_name",
      "seed_eggNOG_ortholog",
      "seed_ortholog_evalue",
      "seed_ortholog_score",
      "Predicted_taxonomic_group",
      "Predicted_protein_name",
      "Gene_Ontology_terms",
      "EC_number",
      "KEGG_ko",
      "KEGG_Pathway",
      "KEGG_Module",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "tax_scope",
      "eggNOG_OGs",
      "bestOG",
      "COG_Functional_Category",
      "eggNOG_free_text_description"
    )
  fread(file,
        col.names = col_names,
        quote = "",
        sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(
      pattern = pattern,
      replacement = "",
      x = basename(file)
    )) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(
      eggNOG_OGs,
      into = c("og", "taxid"),
      sep = "@",
      remove = FALSE
    ) %>%
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
