#!/usr/bin/env Rscript

# check for required packages
list.of.packages <- c("readr", "dplyr","tidyr","stringr","patchwork","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

# load args
args <- commandArgs(trailingOnly=T)
paf_file <- args[1]
refs_comp <- args[2]
sample <- args[3]
gen_frac <- args[4]

# load alignment file
paf <- read_tsv(paf_file, col_names = F) %>%
  rename(TARGET = 6,
         LENGTH = 7,
         START = 8,
         END = 9,
         ALIGN = 11) %>%
  mutate(ASSEMBLY=str_remove_all(TARGET, pattern="(\\d+$)"),
         CONTIG = str_extract(TARGET, "(\\d+$)")) %>%
  select(TARGET, ASSEMBLY, CONTIG, LENGTH, START, END, ALIGN)

# determine missing targets & add back to paf
paf <- read_tsv(refs_comp) %>%
  rename(TARGET = 1,
         LENGTH = 2) %>%
  mutate(ASSEMBLY=str_remove_all(TARGET, pattern="(\\d+$)"),
         CONTIG = str_extract(TARGET, "(\\d+$)"),
         START = 0,
         END = 0,
         ALIGN = 0) %>%
  select(TARGET, ASSEMBLY, CONTIG, LENGTH, START, END, ALIGN) %>%
  filter(!(TARGET %in% paf$TARGET)) %>%
  rbind(paf)

# get alignment depth across each target
expand_align <- function(index){
    subdf <- paf[index,]
    rslt <- data.frame("ASSEMBLY" = subdf$ASSEMBLY, "CONTIG" = subdf$CONTIG, SITE = subdf$START:subdf$END, "LENGTH" = subdf$LENGTH )
    return(rslt)
}
depth <- do.call(rbind, lapply(1:nrow(paf), FUN = expand_align)) %>%
  group_by(ASSEMBLY, CONTIG, SITE, LENGTH) %>%
  count()

# make alignment plots
make_plot <- function(a){
    subdf <- depth %>% filter(ASSEMBLY == a)
    xscales <- subdf %>% select(ASSEMBLY, CONTIG, LENGTH) %>% unique()
    p <- ggplot()+
        geom_histogram(data = subdf, aes(x = SITE, y = n), stat = "identity", color = "black")+
        geom_segment(data = xscales, aes(x = 0, xend = LENGTH, y = 0, yend = 0))+
        facet_grid(ASSEMBLY~CONTIG, scales = "free")+
        theme_bw()+
        theme(strip.text.y = element_text(angle = 0))+
        labs(x = "Site", y = "Depth")
    return(p)
}
plot_list <- paf %>%
  group_by(ASSEMBLY) %>%
  filter(sum(ALIGN) > 0) %>%
  .$ASSEMBLY %>%
  unique()
n_plots <- length(plot_list)
# determine if plots should be made 
if(n_plots > 0){
  plots <- lapply(plot_list, FUN = make_plot)
  p <- wrap_plots(plots, ncol = 1, nrow = n_plots)
  ggsave(plot = p, file = paste0(sample,".ref-genfrac.jpg"), dpi = 300, width = 10, height = 2*n_plots, limitsize = F)
}

# set max depth to 1 & determine genome fraction
ref.summary <- depth %>%
  mutate(n_adj = case_when(n > 0 ~ 1,
                           TRUE ~ 0)) %>%
  group_by(ASSEMBLY, LENGTH) %>%
  summarize(ALIGN = sum(n_adj)-1) %>%
  group_by(ASSEMBLY) %>%
  summarize(GENOME_FRAC = sum(ALIGN) / sum(LENGTH))

# select references
ref.list <- ref.summary %>%
  subset(as.numeric(GENOME_FRAC) > as.numeric(gen_frac)) %>%
  select(ASSEMBLY)
if(nrow(ref.list) == 0){
  ref.list <- data.frame(ASSEMBLY = "none_selected")
}

# write outputs
write.table(x= ref.summary, file = paste0(sample,".ref-summary.csv"), quote = F, sep = ",", row.names = F)
write.table(x= ref.list, file = paste0(sample,".ref-list.csv"), quote = F, sep = ",", row.names = F, col.names = F)
