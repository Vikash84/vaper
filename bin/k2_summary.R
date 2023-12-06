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
k2_output <- args[1]
contig_cov <- args[2]
sample <- args[3]
ncbi_stats <- args[4]

# determine approx. percentages of each taxa observed and output as .tsv
## load contig stats
df.cov <- read.csv(contig_cov, col.names = c("CONTIG", "COVERAGE"))
# load Kraken2 output
df.k2 <- read_tsv(k2_output, col_names = F) %>%
  rename(STATUS = 1,
         CONTIG = 2,
         TAXA_TAXID = 3,
         LENGTH = 4,
         LCA = 5) %>%
  separate(col="TAXA_TAXID", sep = " \\(taxid ", into = c("TAXA", "TAXID")) %>%
  mutate(TAXID = gsub(TAXID, pattern = ")$", replacement = ""))
# load NCBI stats
df.ncbi <- read_tsv(ncbi_stats) %>%
  select(TaxID, `Size (Kb)`) %>%
  rename(TAXID = TaxID, 
         LENGTH = `Size (Kb)`) %>%
  group_by(TAXID) %>%
  summarize(REF_LENGTH = mean(LENGTH))

# combine all datastreams together
df.k2 %>%
  merge(df.cov, by = "CONTIG") %>%
  merge(df.ncbi, by = "TAXID") %>%
  mutate(EST_READS=as.numeric(LENGTH)*as.numeric(COVERAGE)) %>%
  group_by(TAXA, TAXID, REF_LENGTH) %>%
  summarise(GROUP_SUM = sum(EST_READS), TOTAL_LENGTH = sum(LENGTH), MEAN_COVERAGE = mean(COVERAGE)) %>%
  ungroup() %>%
  mutate(TOT_READS = sum(GROUP_SUM),
         EST_PER_ABUND = as.character(100*GROUP_SUM / TOT_READS),
         EST_GEN_FRAC = (TOTAL_LENGTH/1000) / REF_LENGTH) %>%
  select(TAXA, TAXID, EST_PER_ABUND, TOTAL_LENGTH, MEAN_COVERAGE, EST_GEN_FRAC) %>% 
  arrange(desc(EST_PER_ABUND)) %>%
  write.csv(file = paste0(sample,".k2-summary.csv"), quote = F, row.names = F)