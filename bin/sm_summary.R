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
sm_output <- args[1]
sample <- args[2]

# load Sourmash gather output
df.sm <- read.csv(sm_output)
colnames(df.sm)
#  write.csv(file = paste0(sample,".k2-summary.csv"), quote = F, row.names = F)