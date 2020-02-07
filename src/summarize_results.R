#!/bin/env Rscript
#
# Summarizes levels of enrichment for rankings assessed
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

# load fgsea results
dat <- read_parquet(snakemake@input[[1]])

# summarize level of significant enrichment for each covariate
dat %>%
  group_by(field) %>%
  summarize(num_sig = sum(padj < 0.01)) %>%
  write_tsv(snakemake@output[['summary']])

# measure similarity of covariates in terms of functional enrichment
dat_wide <- dat %>%
  pivot_wider(id_cols = pathway, names_from = field, values_from = pval)

# measure correlation of covariates
cor_mat <- cor(dat_wide[, -1], use = 'pairwise.complete.obs')
cor_mat <- as.data.frame(cor_mat)

rownames(cor_mat) <- colnames(cor_mat)

write_tsv(as.data.frame(cor_mat), snakemake@output[['cor_mat']])
