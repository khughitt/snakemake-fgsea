"""
Snakemake functional enrichment pipeline
V. Keith Hughitt
"""
import os
import yaml

# output directory
out_dir = os.path.join(config["output_dir"], config["version"])

# dataset names and filepaths
dataset_names = config["datasets"].keys()

# get a list of gene sets to be analyzed
gmts = os.listdir(config["gene_sets"]["include"])
gene_set_names = [x.replace(".gmt", "") for x in gmts]

rule all:
    input:
        expand(os.path.join(out_dir, "results/summary/{dataset}_summary.tsv"),
               dataset=dataset_names)

rule summarize_results:
    input:
        os.path.join(out_dir, "results/merged/{dataset}.parquet")
    output:
        os.path.join(out_dir, "results/summary/{dataset}_summary.tsv")
    script: "src/summarize_results.R"


rule combine_results:
    input:
        expand(os.path.join(out_dir, "results/{{dataset}}/collections/{gene_set}.parquet"),
               gene_set=gene_set_names)
    output:
        os.path.join(out_dir, "results/merged/{dataset}.parquet")
    script: "src/combine_results.R"

rule run_fgsea:
    input: 
        dataset=os.path.join(out_dir, "input", "{dataset}.feather"),
        gmt=os.path.join(config["gene_sets"]["include"], "{gene_set}.gmt")
    output:
        os.path.join(out_dir, "results/{dataset}/collections/{gene_set}.parquet")
    script: "src/run_fgsea.R"

rule create_dataset_symlinks:
    output:
        os.path.join(out_dir, "input", "{dataset}.feather")
    run:
        os.symlink(config["datasets"][wildcards["dataset"]]["path"], output[0])
