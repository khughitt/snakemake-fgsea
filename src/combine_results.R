#!/bin/env Rscript
#
# Combines enrichment results from multiple collections into a single file
#
suppressMessages(library(arrow))
options(stringsAsFactors = FALSE)

dat <- do.call(rbind, lapply(snakemake@input, read_parquet))

write_parquet(dat, snakemake@output[[1]], compression = 'ZSTD')
