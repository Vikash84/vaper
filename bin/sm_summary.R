#!/usr/bin/env Rscript

# sm_summary.R
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
# load args
args <- commandArgs(trailingOnly=T)
sm_taxa <- args[1]
prefix <- args[2]
# load sourmash data
df <- read_csv(sm_taxa)
# create plot & summary
if(nrow(df) > 0){
    # get fraction of sample that was unclassified
    df.unclass <- df %>%
      filter(rank == "species" & lineage == "unclassified")
    if(nrow(df.unclass) > 0){
      unclass_frac <- df.unclass$fraction
    }else(unclass_frac <- 0)
    # clean up data and calculate percentages
    df <- df %>%
      filter(rank == 'species') %>%
      mutate(classified_fraction = fraction / (1-unclass_frac),
             Species = case_when(classified_fraction >= 0.01 ~ gsub(lineage, pattern = ".*;", replacement = ""),
                                 TRUE ~ "Other")) %>%
      group_by(Species) %>%
      summarize(fraction = sum(fraction), classified_fraction = sum(classified_fraction)) %>%
      ungroup() %>%
      arrange(desc(fraction)) %>%
      mutate(Species = factor(Species, levels = Species))

    # create plot & save
    p <- df %>%
      filter(Species != 'unclassified') %>%
      mutate(ymax = cumsum(100*classified_fraction),
             ymin = c(0, head(ymax, n=-1))) %>%
      ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Species))+
      geom_rect() +
      coord_polar(theta="y") +
      xlim(c(2, 4))+
      theme(axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    ggsave(plot = p, filename = paste0(prefix,".taxa-plot.jpg"), dpi = 300, height = 10, width = 15)

    # create summaryline
    summaryline <- df %>%
      arrange(Species) %>%
      mutate(virus = paste0(round(100*fraction, digits =1),"% ", Species)) %>%
      .$virus %>%
      paste(collapse = "; ") %>%
      gsub(pattern = ",", replacement = "")
    write(x = summaryline, file = paste0(prefix,".taxa-summary.csv"))
}else(write(x = "No Viruses Detected", file = paste0(prefix,".taxa-summary.csv")))