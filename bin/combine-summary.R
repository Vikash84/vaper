#!/usr/bin/env Rscript

version <- "1.0"

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
         "ASSEMBLY_NON_ATCGN",
         "ASSEMBLY_TERMINAL_GAPS",
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
         "SPECIES_SUMMARY")

df.empty <- matrix(ncol = length(col.list), nrow = 1) %>%
  data.frame()
colnames(df.empty) <- col.list

#----- LOAD SUMMARY LINES ----#
# get list of summary lines
files <- list.files("./", pattern = ".csv", full.names = T)
# combine all lines
mergeTables <- function(file){
  read_csv(file) %>%
    mutate_all(as.character) %>%
    return()
}
df <- do.call(bind_rows, lapply(files, FUN=mergeTables)) %>%
  bind_rows(df.empty)


#----- GATHER FINAL METRICS -----#
# depth of coverage & percent read/bp mapping
df <- df %>%
  mutate(ASSEMBLY_EST_DEPTH = round(as.numeric(BASES_MAPPED) / as.numeric(ASSEMBLY_LENGTH), digits = 0),
         PERC_READS_MAPPED = round(100*as.numeric(READS_MAPPED) / as.numeric(TOTAL_READS_CLEAN), digits = 1),
         PERC_BASES_MAPPED = round(100*as.numeric(BASES_MAPPED) / as.numeric(TOTAL_BASES_CLEAN),digits = 1))
# QC status
df <- df %>%
  mutate(ASSEMBLY_QC_REASON = '',
         ASSEMBLY_QC_REASON = case_when(ASSEMBLY_EST_DEPTH < as.numeric(min_depth) ~ 'Low depth of coverage',
                                        TRUE ~ ASSEMBLY_QC_REASON),
         ASSEMBLY_QC_REASON = case_when(ASSEMBLY_GEN_FRAC < as.numeric(min_genfrac) ~ paste(ASSEMBLY_QC_REASON, 'Low genome fraction', sep = "; "),
                                        TRUE ~ ASSEMBLY_QC_REASON),       
         ASSEMBLY_QC = case_when(ASSEMBLY_QC_REASON == '' ~ 'PASS',
                                 TRUE ~ 'FAIL'),
         ASSEMBLY_QC_REASON = case_when(ASSEMBLY_NC == 'mediocre' | ASSEMBLY_NC == 'bad' ~ paste(ASSEMBLY_QC_REASON, "WARNING: Nextclade QC rating was '",ASSEMBLY_NC,"'", sep = '; '),
                                        TRUE ~ ASSEMBLY_QC_REASON),
         ASSEMBLY_QC_REASON = gsub(ASSEMBLY_QC_REASON, pattern = '^; ', replacement = ''),
         ASSEMBLY_QC_REASON = gsub(ASSEMBLY_QC_REASON, pattern = '; $', replacement = ''))
# summarize assembly variants
df <- df %>%
  group_by(ID,SPECIES,SEGMENT) %>%
  mutate(ASSEMBLY_VARIANT = paste0(row_number()," of ", n())) %>%
  ungroup()

# order columns
allCols <- colnames(df)
df <- df %>%
  select(all_of(c(col.list, allCols[!(allCols %in% col.list)]))) %>%
  drop_na(ID)
# save combined summary
write.csv(x=df, file="combined-summary.csv", row.names = F, quote = F)