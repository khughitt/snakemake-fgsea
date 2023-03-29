#
# Measure function enrichment for genes sets from a single gmt file
#
suppressMessages(library(arrow))
suppressMessages(library(GSEABase))
suppressMessages(library(arrow))
suppressMessages(library(fgsea))
suppressMessages(library(tidyverse))

set.seed(1)

# load gene scores
dat <- read_feather(snakemake@input$dataset)

# drop any unneeded columns
dataset_name <- snakemake@wildcards$dataset

dataset_config <- snakemake@config$datasets

if ("exclude" %in% names(dataset_config)) {
  dat <- dat[, !colnames(dat) %in% dataset_config$exclude]
}

# load gene sets
gene_sets <- geneIds(getGmt(snakemake@config$gene_sets$gmt))

fgsea_results <- NULL

# iterate over dataset columns, skipping id column at position 1
for (i in 2:ncol(dat)) {
  # create a named vector of the numeric values to be tested for enrichment
  gene_vals <- deframe(dat[, c(1, i)])

  message(sprintf("Processing: %s...", colnames(dat)[i]))

  # check for problematic edge case where all non-missing p-values are "1"
  if (min(gene_vals, na.rm=TRUE) == 1) {
    stop("Dataset has no significant P-values!")
  }

  # transform values, if requested
  if ("transform" %in% names(dataset_config)) {
    if (dataset_config[["transform"]] == "-log10") {
      gene_vals <- -log10(pmax(gene_vals, 1E-20))
    } else {
      stop("Invalid dataset transform specified!")
    }
  }

  # drop missing values
  gene_vals <- gene_vals[!is.na(gene_vals)]

  cname <- colnames(dat)[i]

  set.seed(1)

  # run fgsea on the gene scores
  res <- fgsea(gene_sets, 
               stats = gene_vals,
               scoreType = snakemake@config$fgsea$score_type,
               eps = as.numeric(snakemake@config$fgsea$eps),
               minSize = snakemake@config$gene_sets$min_size,
               maxSize = snakemake@config$gene_sets$max_size,
               nproc = snakemake@config$num_threads) %>%
    select(-leadingEdge) %>%
    arrange(pval)

  message("Done!")

  entry <- cbind(dataset = dataset_name, field = cname,
                 gene_set = snakemake@config$gene_sets$gmt, res)
  fgsea_results <- rbind(fgsea_results, entry)
}

# drop entries for which no p-value could be computed
fgsea_results <- fgsea_results[complete.cases(fgsea_results), ]

write_feather(fgsea_results, snakemake@output[[1]])
