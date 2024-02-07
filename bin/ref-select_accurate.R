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
paf_file <- args[1]
refs_comp <- args[2]
sample <- args[3]
gen_frac <- args[4]


# load alignment file
paf <- read_tsv(paf_file, col_names = F) %>%
  rename(CONTIG = 6,
         LENGTH = 7,
         ALIGN = 11) %>%
  mutate(ASSEMBLY=gsub(CONTIG, pattern="[0-9]$",replacement=""))

# load reference contig lengths
contig_lengths <- read_tsv(refs_comp) %>%
  rename(CONTIG = 1,
         LENGTH = 2) %>%
  mutate(ASSEMBLY=gsub(CONTIG, pattern="[0-9]$",replacement="")) %>%
  group_by(ASSEMBLY) %>%
  summarise(TOT_LENGTH = sum(LENGTH))

# sum alignments
align.summary <- paf %>%
  group_by(ASSEMBLY) %>%
  summarise(TOT_ALIGN = sum(ALIGN)) %>%
  merge(contig_lengths, by = "ASSEMBLY") %>%
  mutate(GENOME_FRAC = TOT_ALIGN / TOT_LENGTH)

# select references
ref.list <- align.summary %>%
  subset(GENOME_FRAC > gen_frac) %>%
  select(ASSEMBLY)
if(nrow(ref.list) == 0){
  ref.list <- data.frame(ASSEMBLY = "none_selected")
}

# write outputs
write.table(x= align.summary, file = paste0(sample,".ref-summary.csv"), quote = F, sep = ",", row.names = F)
write.table(x= ref.list, file = paste0(sample,".ref-list.csv"), quote = F, sep = ",", row.names = F, col.names = F)
