#!/usr/bin/env Rscript

version <- "1.1"

# combine-summary.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- LIBRARIES -----#
# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)

#----- ARGUMENTS -----#
# load args
args <- commandArgs(trailingOnly=T)
min_depth <- args[1]
min_genfrac <- args[2]

#---- VERSION ----#
if(args[1] == "version"){
  cat(version, sep = "\n")
  quit(status=0)
}

#----- FUNCTIONS -----#
qcCheck <- function(ref, gf, gf_threshold, depth, depth_threshold){
  if(is.na(gf) || is.na(depth)){
    if(!is.na(ref)){
      status <- 'FAIL'
      reason <- 'No assembly'
    }else{
      # No assembly
      status <- NA_character_
      reason <- NA_character_
    }
  }else{
    # Default status
    status <- 'PASS'
    reason_gf <- NA_character_
    reason_depth <- NA_character_

    # Genome fraction
    if(as.numeric(gf) < as.numeric(gf_threshold)){
      status <- 'FAIL'
      reason_gf <- paste0('Genome fraction below ',gf_threshold)
    }
    # Depth of coverage
    if(as.numeric(depth) < as.numeric(depth_threshold)){
      status <- 'FAIL'
      reason_depth <- paste0('Depth of coverage below ',depth_threshold,"X")
    }
    reason <- c(reason_gf, reason_depth) %>%
      na.omit() %>%
      unique() %>%
      paste(collapse = ';')
    if(reason == ''){ reason <- NA_character_ }
  }
  
  return(list(status, reason))
}

#----- CREATE EMPTY DATAFRAME -----#
# set expected columns & create empty dataframe
col.list <- c("ID",
         "REFERENCE",
         "SPECIES",
         "SEGMENT",
         "ASSEMBLY_VARIANT",
         "ASSEMBLY_QC",
         "ASSEMBLY_QC_REASON",
         "ASSEMBLY_LENGTH",
         "REF_LENGTH",
         "ASSEMBLY_EST_DEPTH",
         "ASSEMBLY_GEN_FRAC",
         "ASSEMBLY_SUBS",
         "ASSEMBLY_DEL",
         "ASSEMBLY_INS",
         "ASSEMBLY_MISSING",
         "ASSEMBLY_NON_ACGTN",
         "ASSEMBLY_TERMINI_GAPS",
         "PERC_READS_MAPPED",
         "PERC_BASES_MAPPED",
         "READS_MAPPED",
         "BASES_MAPPED",
         "MEAN_MAPPED_READ_LENGTH",
         "MEAN_MAPPED_READ_QUALITY",
         "TOTAL_READS_CLEAN",
         "TOTAL_BASES_CLEAN",
         "READ1_MEAN_LENGTH_CLEAN",
         "READ2_MEAN_LENGTH_CLEAN",
         "Q30_RATE_CLEAN",
         "TOTAL_READS_RAW",
         "TOTAL_BASES_RAW",
         "READ1_MEAN_LENGTH_RAW",
         "READ2_MEAN_LENGTH_RAW",
         "Q30_RATE_RAW",
         "SPECIES_SUMMARY",
         "REF_INFO")

df.empty <- matrix(ncol = length(col.list), nrow = 1) %>%
  data.frame()
colnames(df.empty) <- col.list

#----- LOAD SUMMARY LINES ----#
# get list of summary lines
files <- list.files("./", pattern = ".csv", full.names = T)
# combine all lines
mergeTables <- function(file){
  read_csv(file) %>%
    mutate(across(where(is.character), ~na_if(., "null"))) %>%
    mutate_all(as.character) %>%
    return()
}
df <- do.call(bind_rows, lapply(files, FUN=mergeTables)) %>%
  bind_rows(df.empty) %>%
  drop_na(ID)
#----- GATHER FINAL METRICS -----#
# depth of coverage & percent read/bp mapping
df <- df %>%
  mutate(ASSEMBLY_EST_DEPTH = round(as.numeric(BASES_MAPPED) / as.numeric(ASSEMBLY_LENGTH), digits = 0),
         PERC_READS_MAPPED = round(100*as.numeric(READS_MAPPED) / as.numeric(TOTAL_READS_CLEAN), digits = 1),
         PERC_BASES_MAPPED = round(100*as.numeric(BASES_MAPPED) / as.numeric(TOTAL_BASES_CLEAN),digits = 1))
# QC status
df <- df %>%
  group_by(ID,REFERENCE) %>%
  mutate(ASSEMBLY_QC = list(qcCheck(REFERENCE,ASSEMBLY_GEN_FRAC,min_genfrac,ASSEMBLY_EST_DEPTH,min_depth)),
         ASSEMBLY_QC_REASON = unlist(ASSEMBLY_QC)[2],
         ASSEMBLY_QC = unlist(ASSEMBLY_QC)[1]) %>%
  ungroup() %>%
  select(-any_of(c('ASSEMBLY_NC')))
# summarize assembly variants
df <- df %>%
  group_by(ID,SPECIES,SEGMENT) %>%
  mutate(ASSEMBLY_VARIANT = case_when( !is.na(SPECIES) ~ paste0(row_number()," of ", n()),
                                        TRUE ~ NA_character_ )) %>%
  ungroup()

# order columns
allCols <- colnames(df)
df <- df %>%
  select(all_of(c(col.list, allCols[!(allCols %in% col.list)]))) %>%
  select_if(~ !all(is.na(.))) %>%
  drop_na(ID)
# save combined summary
write.csv(x=df, file="combined-summary.csv", row.names = F, quote = F)