#
# Measure function enrichment for genes sets from a single gmt file
#
suppressMessages(library(arrow))
suppressMessages(library(GSEABase))
suppressMessages(library(arrow))
suppressMessages(library(fgsea))
suppressMessages(library(tidyverse))

set.seed(1)

# configuration parameters
min_size <- snakemake@config$gene_sets$min_size
max_size <- snakemake@config$gene_sets$max_size
fgsea_nperm <- snakemake@config$gene_sets$fgsea_nperm

# load gene scores
dat <- read_feather(snakemake@input$dataset)

# drop any unneeded columns
dataset_name <- snakemake@wildcards$dataset

dataset_config <- snakemake@config$datasets

if ("exclude" %in% names(dataset_config)) {
  dat <- dat[, !colnames(dat) %in% dataset_config$exclude]
}

# load gene sets
gene_sets = geneIds(getGmt(snakemake@input$gmt))
gene_sets <- gene_sets[gene_sets != '']

# remove any leading or trailing whitespace in gene set names;
# encountered at least one instance of this ("HEK 293 T-rex" in NCI-60 gene set
# gene_sets)
names(gene_sets) <- trimws(names(gene_sets))

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
gene_sets <- lapply(gene_sets, function(x) {
  sub(',\\d+\\.\\d+$', '', x)
})

# exclude any gene sets which are either too large or too small
set_sizes <- lapply(gene_sets, length)
mask <- (set_sizes >= min_size) & (set_sizes <= max_size)
gene_sets <- gene_sets[mask]

# work-around for possible corrupted entry in DSigDB GMT (line 6027);
# issue has been reported and manually corrected in the 1.0 currently in use, but
# leaving this code in for now in case it accidentally gets reverted.
if (snakemake@wildcards[['gene_set']] == "DSigDB_All") {
  ind <- which(names(gene_sets) == '')
  gene_sets <- gene_sets[-ind]
}

fgsea_results <- NULL

# iterate over dataset columns, skipping id column at position 1
for (i in 2:ncol(dat)) {
  # create a named vector of the numeric values to be tested for enrichment
  gene_vals <- deframe(dat[, c(1, i)])

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
  res <- fgsea(gene_sets, stats = gene_vals, nperm = fgsea_nperm, 
               nproc = snakemake@config$num_threads) %>%
    select(-leadingEdge) %>%
    arrange(pval)

  entry <- cbind(dataset = dataset_name, field = cname, 
                 gene_set = snakemake@wildcards$gene_set, res)
  fgsea_results <- rbind(fgsea_results, entry)
}

# drop entries for which no p-value could be computed
fgsea_results <- fgsea_results[complete.cases(fgsea_results), ]

# write_parquet(fgsea_results, snakemake@output[[1]], compression = 'ZSTD')
write_feather(fgsea_results, snakemake@output[[1]])
