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
    mutate(ASSEMBLY_EST_DEPTH = BASES_MAPPED / ASSEMBLY_LENGTH)
}
# save combined summary
write.csv(x=df, file="combined-summary.csv", quote = F, row.names = F)