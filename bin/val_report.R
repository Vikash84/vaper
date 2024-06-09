#!/usr/bin/env Rscript

# val_report.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)

# calculate global metrics
global_metrics <- function(results_file, pairs_file, metric){
    if((file.exists(results_file) && file.exists(pairs_file))){
        df <- read_csv(results_file)
        results <- df[,ncol(df)] %>% unlist()
        n_missing_extra <- read_tsv(pairs_file, col_names = F) %>%
          rename(seq1=1,
                 seq2=2) %>%
          filter(seq1 == "null" | seq2 == "null") %>%
          nrow()
        results <- c(results, rep(0, n_missing_extra))
        n <- length(results)
        p_99 <- 100*sum(results >= 99)/n
        p_95 <- 100*sum(results >= 95)/n
        p_90 <- 100*sum(results >= 90)/n
        return(data.frame("metric" = metric, "min"=min(results), "max" = max(results), "mean" = mean(results), "sdev" = sd(results), "med" = median(results), "n" = n, "n_error" = n_missing_extra, "p_99" = p_99, "p_95" = p_95, "p_90"=p_90))
    }else(return(data.frame("metric" = metric)))
}

acc <- global_metrics("accuracy_results.csv","accuracy_pairs.csv", "Accuracy")
inter <- global_metrics("precision_inter_results.csv","precision_inter_pairs.csv", "Inter-Assay Reproducility")
intra <- global_metrics("precision_intra_results.csv","precision_intra_pairs.csv", "Intra-Assay Reproducility")

final <- do.call(bind_rows, list(acc, inter, intra)) %>%
  mutate()
write.csv(x=final, file = "validation-report.csv", quote = F, row.names = F)