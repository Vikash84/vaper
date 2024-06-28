#!/usr/bin/env Rscript

# val_report.R
# Author: Jared Johnson, jared.johnson@doh.wa.gov

# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# calculate global metrics
global_metrics <- function(results_file, pairs_file, metric){
    if((file.exists(results_file) && file.exists(pairs_file))){
        df <- read_csv(results_file)
        results <- df[,ncol(df)] %>% unlist()
        df.pairs <- read_tsv(pairs_file, col_names = F) %>%
          rename(seq1=1,
                 seq2=2)
        n_missing <- df.pairs %>%        
          filter(seq1 == "null") %>%
          nrow()
        n_extra <- df.pairs %>%        
          filter(seq2 == "null") %>%
          nrow()
        n <- length(results)
        p_99 <- 100*sum(results >= 99)/n
        p_95 <- 100*sum(results >= 95)/n
        p_90 <- 100*sum(results >= 90)/n
        p <- df[,ncol(df)] %>%
          rename(value = 1) %>%
          ggplot(aes(x=value))+
            geom_histogram()+
            theme_bw()
        ggsave(plot = p, filename = paste0(metric,".jpg"),dpi=300, width = 5, height = 3)
        return(data.frame("metric" = metric, "min"=min(results), "max" = max(results), "mean" = mean(results), "sdev" = sd(results), "med" = median(results), "n" = n, "n_missing" = n_missing, "n_extra" = n_extra, "p_99" = p_99, "p_95" = p_95, "p_90"=p_90))
    }else(return(data.frame("metric" = metric)))
}

acc <- global_metrics("accuracy_results.csv","accuracy_pairs.tsv", "Accuracy")
inter <- global_metrics("precision_inter_results.csv","precision_inter_pairs.tsv", "Inter-Assay Reproducility")
intra <- global_metrics("precision_intra_results.csv","precision_intra_pairs.tsv", "Intra-Assay Reproducility")

final <- do.call(bind_rows, list(acc, inter, intra)) %>%
  mutate()
write.csv(x=final, file = "validation-report.csv", quote = F, row.names = F)