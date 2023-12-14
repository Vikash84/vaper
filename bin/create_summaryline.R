#!/usr/bin/env Rscript

library(tidyverse)

# load args
args <- commandArgs(trailingOnly=T)
fastp2tbl <- args[1]
k2_summary <- args[2]
samtoolstats2tbl <- args[3]
assembly_stats <- args[4]
sample <- args[5]
ref <- args[6]
snvs <- args[7]

#----- Sample ID & Reference
df.summaryline <- data.frame(ID = sample, REFERENCE = ref)

#----- Full Read Stats -----#
df.fastp2tbl <- read_csv(fastp2tbl)
df.summaryline <- cbind(df.summaryline, df.fastp2tbl)

#----- Kraken2 Summary -----#
# load summary
df.k2 <- read.csv(k2_summary) %>%
  arrange(desc(EST_PER_ABUND))
# determine main taxa - i.e., those present at >= 1% relative abundance
main_taxa <- df.k2 %>%
  filter(EST_PER_ABUND >= 1) %>%
  .$TAXA %>%
  paste(collapse = "; ")
# determine other taxa - i.e., those preset at < 1% relative abundance
other_taxa <- df.k2 %>%
  filter(EST_PER_ABUND < 1) %>%
  .$TAXA %>%
  paste(collapse = "; ")
# combine into a table
df.k2_main_other <- data.frame(MAIN_TAXA_KRAKEN2 = main_taxa, OTHER_TAXA_KRAKEN2 = other_taxa)
df.summaryline <- cbind(df.summaryline, df.k2_main_other)

#----- Mapped Reads Stats -----#
if(file.exists(samtoolstats2tbl)){
  df.samtoolstats2tbl <- read.csv(samtoolstats2tbl)
  df.summaryline <- cbind(df.summaryline, df.samtoolstats2tbl)
}else(cat("\nNo mapping stats provided\n"))

#----- Assembly Stats -----#
if(file.exists(assembly_stats)){
  cat(snvs)
  df.assembly_stats <- read.csv(assembly_stats) %>%
    mutate(assembly_snvs = snvs,
           assembly_nuc_id = (assembly_length - assembly_snvs) / assembly_length)
  df.summaryline <- cbind(df.summaryline, df.assembly_stats)
}else(cat("\nNo assembly stats provided\n"))

#----- Create Summaryline -----#
df.summaryline <- df.summaryline %>% 
      rename_with(toupper)
write.csv(x=df.summaryline, file = "summaryline.csv", quote = F, row.names = F)
