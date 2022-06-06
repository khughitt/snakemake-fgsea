#!/bin/env Rscript
#
# Generates a Jaccard index-based gene set similarity matrix from a collection of gmt
# files
#
library(arrow)
library(GSEABase)

# load gene sets from gmt files
infiles <- Sys.glob(file.path(snakemake@config$gene_sets$include, "*.gmt"))

gene_sets <- lapply(infiles, function(x) {
  gs <- geneIds(getGmt(x))

  # remove any gene sets with missing names
  gs <- gs[gs != '']

  # remove any leading or trailing whitespace in gene set names
  names(gs) <- trimws(names(gs))

  # remove gene set weights, if present (e.g. "ANXA1,1.0" -> "ANXA1")
  gs <- lapply(gs, function(x) {
    sub(',\\d+\\.\\d+$', '', x)
  })

  gs
})
names(gene_sets) <- sub('.gmt', '', basename(infiles))

# add collection name as a prefix to gene set identifiers
for (collection in names(gene_sets)) {
  names(gene_sets[[collection]]) <- paste0(collection, '_', names(gene_sets[[collection]]))
}

# create a mapping from the combined <collection>_<geneset> identifiers to their
# corresponding collection and gene set identifiers
mapping <- NULL

for (collection in names(gene_sets)) {
  combined_ids <- paste0(collection, '_', names(gene_sets[[collection]]))

  mapping <- rbind(mapping, cbind(combined_id = combined_ids, collection, 
                                  gene_set = names(gene_sets[[collection]])))
}
mapping <- as.data.frame(mapping)
write_feather(mapping, snakemake@output$mapping)

# collapse entries from multiple files into a single list 
gene_sets <- unlist(gene_sets, recursive = FALSE)

# generate a matrix with to store gene set overlap information
sim_mat <- matrix(NA, nrow = length(gene_sets), ncol = length(gene_sets))

# iterate over gene set pairs and compute overlap (jaccard indices)
for (i in 1:length(gene_sets)) {
  genes1 <- gene_sets[[i]]

  for (j in 1:length(gene_sets)) {
    # self-comparisons
    if (i == j) {
      sim_mat[i, j] <- 1
      next
    }

    # skip symmetric values already computed
    if (!is.na(sim_mat[i, j])) {
      next
    }

    # otherwise, compute overlap
    genes2 <- gene_sets[[j]]

    sim_mat[i, j] <- length(intersect(genes1, genes2)) / length(union(genes1, genes2))
    sim_mat[j, i] <- sim_mat[i, j]
  }
}
rownames(sim_mat) <- names(gene_sets)
colnames(sim_mat) <- names(gene_sets)

# store similarity matrix
sim_mat <- as.data.frame(sim_mat)

write_feather(sim_mat, snakemake@output$sim_mat)
