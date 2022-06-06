#
# Generates an undirected network of gene sets with node weights corresponding 
# to fgsea enrichment significance and edge weights corresponding to the similarity of 
# the two connected gene sets (based on correlation of overlapping genes with
# all gene sets)
#
suppressMessages(library(arrow))
suppressMessages(library(igraph))
suppressMessages(library(tidyverse))

# load fgsea results
fgsea_res <- read_feather(snakemake@input[[1]]) %>%
  filter(field == snakemake@config$network$target_field)

# load gene set correlation matrix
cor_mat <- readRDS(snakemake@config$network$correlations)

# load gene set clustering
gene_set_clusters <- read_tsv(snakemake@config$network$clusters)

# remove gene sets with low significance
fgsea_res <- fgsea_res %>%
  filter(padj < as.numeric(snakemake@config$network$padj_cutoff))

mask <- rownames(cor_mat) %in% fgsea_res$pathway
cor_mat <- cor_mat[mask, mask]

gene_set_clusters <- gene_set_clusters[mask, ]

# reorder fgsea results to match correlation matrix ordering
ind <- match(rownames(cor_mat), fgsea_res$pathway)
fgsea_res <- fgsea_res[ind, ]

# remove edges corresponding to gene sets with little similarity
mask <- cor_mat <= as.numeric(snakemake@config$network$min_cor)
cor_mat[mask] <- 0

# create igraph network
g <- graph.adjacency(cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE)

# add vertex properties
g <- set.vertex.attribute(g, 'collection', value = fgsea_res$gene_set)
g <- set.vertex.attribute(g, 'cluster', value = gene_set_clusters$cluster)
g <- set.vertex.attribute(g, 'pval', value = fgsea_res$pval)
g <- set.vertex.attribute(g, 'padj', value = fgsea_res$padj)
g <- set.vertex.attribute(g, 'size', value = fgsea_res$size)

# Save graph to a file
write.graph(g, snakemake@output[[1]], format = "graphml")

