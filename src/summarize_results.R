#!/bin/env Rscript
#
# Summarizes levels of enrichment for rankings assessed
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

dat <- read_parquet(snakemake@input[[1]])

dat %>%
  group_by(field) %>%
  summarize(num_sig = sum(padj < 0.01)) %>%
  write_tsv(snakemake@output[[1]])
