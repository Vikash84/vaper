#!/usr/bin/env Rscript

# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)

# get list of summary lines
files <- list.files("./", pattern = ".csv", full.names = T)
print(files)
# combine all lines
df <- do.call(bind_rows, lapply(files, FUN=read.csv))
# calculate depth of coverage, if any of the samples mapped
if("BASES_MAPPED" %in% colnames(df)){
  df <- df %>%
    mutate(ASSEMBLY_EST_DEPTH = round(BASES_MAPPED / ASSEMBLY_LENGTH, digits = 0))
}
# order columns
df <- df %>%
  select(ID,
         REFERENCE,
         ASSEMBLY_QC,
         ASSEMBLY_LENGTH,
         ASSEMBLY_EST_DEPTH,
         ASSEMBLY_GEN_FRAC,
         ASSEMBLY_SUBS,
         ASSEMBLY_DEL,
         ASSEMBLY_INS,
         ASSEMBLY_FRAMESHIFTS,
         ASSEMBLY_NON_ATCG,
         READS_MAPPED,
         BASES_MAPPED,
         MEAN_MAPPED_READ_LENGTH,
         MEAN_MAPPED_READ_QUALITY,
         PERC_PAIRED_MAPPED_READS,
         TOTAL_READS_CLEAN,
         TOTAL_BASES_CLEAN,
         READ1_MEAN_LENGTH_CLEAN,
         READ2_MEAN_LENGTH_CLEAN,
         Q30_RATE_CLEAN,
         TOTAL_READS_RAW,
         TOTAL_BASES_RAW,
         READ1_MEAN_LENGTH_RAW,
         READ2_MEAN_LENGTH_RAW,
         Q30_RATE_RAW,
         SPECIES_SUMMARY
         )
# save combined summary
write.csv(x=df, file="combined-summary.csv", quote = F, row.names = F)