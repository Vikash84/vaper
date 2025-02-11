#!/usr/bin/env Rscript

version <- "1.1"

# summaryline.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)

# load args
args <- commandArgs(trailingOnly=T)
fastp2tbl <- args[1]
sm_summary <- args[2]
samtoolstats2tbl <- args[3]
nextclade <- args[4]
sample <- args[5]
ref <- args[6]
refsheet <- args[7]

#---- VERSION ----#
if(args[1] == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#----- Functions -----#
getRefName <- function(filename){
  tokens <- filename %>%
    basename() %>%
    strsplit('\\.') %>% 
    unlist()
  n_tokens <- length(tokens)
  index <- n_tokens - 1
  if(tokens[n_tokens] == 'gz'){
    index <- n_tokens - 2
  }
  name <- paste(tokens[1:index], collapse = '.')
  return(name)
}
#----- Load Reference Info -----#
df.refs <- read_csv(refsheet) %>%
  rename_all(toupper) %>%
  group_by(ASSEMBLY) %>%
  mutate(REFERENCE = getRefName(ASSEMBLY)) %>%
  ungroup() %>%
  select(-ASSEMBLY) %>%
  drop_na(REFERENCE) %>%
  unique()
    
#----- Sample ID & Reference
df.summaryline <- data.frame(ID = sample, REFERENCE = ref) %>%
  left_join(df.refs, by = "REFERENCE") %>%
  mutate(REFERENCE = case_when(REFERENCE == 'No_Reference' ~ NA_character_,
                               TRUE ~ REFERENCE))

#----- Full Read Stats -----#
df.fastp2tbl <- read_csv(fastp2tbl)
df.summaryline <- cbind(df.summaryline, df.fastp2tbl)

#----- Mapped Reads Stats -----#
if(file.exists(samtoolstats2tbl)){
  df.samtoolstats2tbl <- read.csv(samtoolstats2tbl)
  df.summaryline <- cbind(df.summaryline, df.samtoolstats2tbl)
}else(cat("\nNo mapping stats provided\n"))

#----- Nextclade -----#
if(file.exists(nextclade)){
  df.nextclade <- read_csv(nextclade, col_names = c("ASSEMBLY_LENGTH", 
                                                    "REF_LENGTH", 
                                                    "ASSEMBLY_SUBS",
                                                    "ASSEMBLY_DEL", 
                                                    "ASSEMBLY_INS", 
                                                    "ASSEMBLY_MISSING", 
                                                    "ASSEMBLY_TERMINI_GAPS", 
                                                    "ASSEMBLY_NON_ACGTN", 
                                                    "ASSEMBLY_GEN_FRAC",
                                                    "ASSEMBLY_NC" ))
  df.summaryline <- cbind(df.summaryline, df.nextclade)
}else(cat("\nNo assembly stats provided\n"))

#---- Sourmash Species Summary ----#
df.summaryline <- df.summaryline %>%
  mutate(SPECIES_SUMMARY = read_csv(sm_summary, col_names = F)$X1)

#----- Create Summaryline -----#
# make ID always show up first
column.names <- df.summaryline %>% 
  select(-ID) %>%
  colnames()
df.summaryline <- df.summaryline[,c("ID", column.names)] %>%
  rename_with(toupper)
write.csv(x=df.summaryline, file = "summaryline.csv", quote = F, row.names = F)
