#!/usr/bin/env Rscript

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

#----- Sample ID & Reference
df.summaryline <- data.frame(ID = sample, REFERENCE = ref)

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
  df.nextclade <- read_tsv(nextclade) %>%
      select(qc.overallStatus,
             totalSubstitutions,
             totalDeletions,
             totalInsertions,
             totalFrameShifts,
             totalNonACGTNs,
             coverage,
             ASSEMBLY_LENGTH) %>%
      rename(ASSEMBLY_QC = qc.overallStatus,
             ASSEMBLY_SUBS = totalSubstitutions,
             ASSEMBLY_DEL = totalDeletions,
             ASSEMBLY_INS = totalInsertions,
             ASSEMBLY_FRAMESHIFTS = totalFrameShifts,
             ASSEMBLY_NON_ATCG = totalNonACGTNs,
             ASSEMBLY_GEN_FRAC = coverage)
  df.summaryline <- cbind(df.summaryline, df.nextclade)
}else(cat("\nNo assembly stats provided\n"))

#---- Sourmash Species Summary ----#
df.summaryline <- df.summaryline %>%
  mutate(SPECIES_SUMMARY = gsub(readLines(sm_summary), pattern = ";$", replacement = ""))

#----- Create Summaryline -----#
# make ID always show up first
column.names <- df.summaryline %>% 
  select(-ID) %>%
  colnames()
df.summaryline <- df.summaryline[,c("ID", column.names)] %>%
  rename_with(toupper)
write.csv(x=df.summaryline, file = "summaryline.csv", quote = F, row.names = F)
