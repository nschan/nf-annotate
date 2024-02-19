#!/usr/bin/env Rscript

# This script reimplements what IPS2FpGs.sh does:
# It filters tables from Pfam and Superfamily to identify R genes

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))

args <- commandArgs(trailingOnly = TRUE)

# Args are:
# file1 (e.g. Pfam)
# file2 (e.g. Superfam)
# more files
# outfile (Basename only)

# Define domain - pfam or superfam relationships
## Currently, these use what is defined by HRP

domains_interpro <-
  data.frame(
    Class = c(
      "C",
      "T",
      "T",
      "R",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "N"
    ),
    Family = c(
      "IPR041118",
      "IPR000157",
      "IPR035897",
      "IPR044974",
      "IPR001611",
      "IPR003591",
      "IPR006553",
      "IPR011713",
      "IPR013101",
      "IPR013210",
      "IPR032675",
      "IPR002182"
    )
  )

domains_pfam <-
  data.frame(
    Class = c(
      "C",
      "T",
      "T",
      "T",
      "R",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "N",
      "N"
    ),
    Family = c(
      "PF18052",
      "PF01582",
      "PF13676",
      "SSF52200",
      "PF05659",
      "SSF52047",
      "PF00560",
      "PF08263",
      "PF12799",
      "PF13516",
      "PF13855",
      "PF14580",
      "SSF52058",
      "SSF52540",
      "PF00931"
    )
  )

# Subset into individual classes

TIR_interpro <-  filter(domains_interpro, Class == "T")
COIL_interpro <- filter(domains_interpro, Class == "C")
RPW8_interpro <- filter(domains_interpro, Class == "R")
LRR_interpro <-  filter(domains_interpro, Class == "L")
NB_interpro <-  filter(domains_interpro, Class == "N")
TIR_pfam <-  filter(domains_pfam, Class == "T")
COIL_pfam <- filter(domains_pfam, Class == "C")
RPW8_pfam <- filter(domains_pfam, Class == "R")
LRR_pfam <-  filter(domains_pfam, Class == "L")
NB_pfam <-  filter(domains_pfam, Class == "N")

# Read in files

lapply(args[1:(length(args)-1)],
       \(x) read_tsv(x, col_names = F, show_col_types = FALSE)) %>%
  bind_rows() %>%
  magrittr::set_colnames(
    c(
      "ID",
      "hash",
      "Length",
      "Analysis",
      "Family",
      "Description",
      "Start",
      "End",
      "Score",
      "Status",
      "Date",
      "InterPro_Accession",
      "InterPro_description",
      "extra1",
      "extra2"
    )
  ) %>% 
  # Check TIR, COIL, RPW, NB and LRR presence
  mutate(has_TIR = case_when(InterPro_Accession %in% TIR_interpro$Family ~ TRUE,
                             Family %in% TIR_pfam$Family ~ TRUE,
                             TRUE ~ FALSE),
         has_TIR_evidence = case_when(InterPro_Accession %in% TIR_pfam$Family ~ InterPro_Accession,
                                      Family %in% TIR_pfam$Family ~ Family,
                                      TRUE ~ NA_character_),
         has_COIL = case_when(InterPro_Accession %in% COIL_interpro$Family ~ TRUE,
                              Family %in% COIL_pfam$Family ~ TRUE,
                              TRUE ~ FALSE),
         has_COIL_evidence = case_when(InterPro_Accession %in% COIL_interpro$Family ~ InterPro_Accession,
                                       Family %in% COIL_pfam$Family ~ Family,
                                       TRUE ~ NA_character_),
         has_RPW8 = case_when(InterPro_Accession %in% RPW8_interpro$Family ~ TRUE,
                              Family %in% RPW8_pfam$Family ~ TRUE,
                              TRUE ~ FALSE),
         has_RPW8_evidence = case_when(InterPro_Accession %in% RPW8_interpro$Family ~ InterPro_Accession,
                                       Family %in% RPW8_pfam$Family ~ Family,
                                       TRUE ~ NA_character_),
         has_NB = case_when(InterPro_Accession %in% NB_interpro$Family ~ TRUE,
                            Family %in% NB_pfam$Family ~ TRUE,
                            TRUE ~ FALSE),
         has_NB_evidence = case_when(InterPro_Accession %in% NB_interpro$Family ~ InterPro_Accession,
                                     Family %in% NB_pfam$Family ~ Family,
                                     TRUE ~ NA_character_),
         has_LRR = case_when(InterPro_Accession %in% LRR_interpro$Family ~ TRUE,
                             Family %in% LRR_pfam$Family ~ TRUE,
                             TRUE ~ FALSE),
         has_LRR_evidence = case_when(InterPro_Accession %in% LRR_interpro$Family ~ InterPro_Accession,
                                      Family %in% LRR_pfam$Family ~ Family,
                                      TRUE ~ NA_character_)) %>% 
  # Make uniques
  group_by(ID) %>% 
  summarize(has_TIR = case_when(sum(has_TIR)  > 0 ~ TRUE, TRUE ~ FALSE),
            has_COIL = case_when(sum(has_COIL) > 0 ~ TRUE, TRUE ~ FALSE),
            has_RPW8 = case_when(sum(has_RPW8) > 0 ~ TRUE, TRUE ~ FALSE),
            has_NB = case_when(sum(has_NB) > 0 ~ TRUE, TRUE ~ FALSE),
            has_LRR = case_when(sum(has_LRR) > 0 ~ TRUE, TRUE ~ FALSE),
            TIR_evidence = str_flatten(has_TIR_evidence %>% {.[!is.na(.)]}, collapse = " "),
            COIL_evidence = str_flatten(has_COIL_evidence %>% {.[!is.na(.)]}, collapse = " "),
            RPW8_evidence = str_flatten(has_RPW8_evidence %>% {.[!is.na(.)]}, collapse = " "),
            NB_evidence = str_flatten(has_NB_evidence %>% {.[!is.na(.)]}, collapse = " "),
            LRR_evidence = str_flatten(has_LRR_evidence %>% {.[!is.na(.)]}, collapse = " ")
            ) %>% 
  # Filter of NB-LRR
  filter(has_RPW8 | (has_NB & has_LRR )) %>% 
  # Annotate
  mutate(Type = case_when(
    has_TIR ~ "TNL",
    has_COIL ~ "CNL",
    has_RPW8 ~ "NLR",
    TRUE ~ "NL"
  )) %>% 
  ungroup() %T>%
  # Write larger table
  write_tsv(paste0(args[length(args)],"_NLR_table.tsv")) %>% 
  dplyr::select(ID, Type) %>% 
  # Write relevant table
  write_tsv(paste0(args[length(args)],"_NLR_genes.tsv"))
         