"""
Snakemake functional enrichment pipeline
V. Keith Hughitt
"""
import os
import glob
import yaml
from pathlib import Path

# output directory
out_dir = os.path.join(config["output_dir"], config["version"])

# dataset names and filepaths
dset_paths = []

for input_path in config["datasets"]["paths"]:
    dset_paths = dset_paths + glob.glob(input_path)

dset_names = [Path(x).stem for x in dset_paths]

dset_mapping = {dset_names[i]:dset_paths[i] for i in range(len(dset_names))}

# get a list of gene sets to be analyzed
gmts = os.listdir(config["gene_sets"]["include"])
gene_set_names = [x.replace(".gmt", "") for x in gmts]

rule all:
    input:
        expand(os.path.join(out_dir, "results/summary/{dataset}_cor_mat.tsv"),
               dataset=dset_names)

rule summarize_results:
    input:
        os.path.join(out_dir, "results/merged/{dataset}.feather")
    output:
        summary=os.path.join(out_dir, "results/summary/{dataset}_summary.tsv"),
        cor_mat=os.path.join(out_dir, "results/summary/{dataset}_cor_mat.tsv")
    script: "src/summarize_results.R"


rule combine_results:
    input:
        expand(os.path.join(out_dir, "results/{{dataset}}/collections/{gene_set}.feather"),
               gene_set=gene_set_names)
    output:
        os.path.join(out_dir, "results/merged/{dataset}.feather")
    script: "src/combine_results.R"

rule run_fgsea:
    input: 
        dataset=os.path.join(out_dir, "input", "{dataset}.feather"),
        gmt=os.path.join(config["gene_sets"]["include"], "{gene_set}.gmt")
    output:
        os.path.join(out_dir, "results/{dataset}/collections/{gene_set}.feather")
    script: "src/run_fgsea.R"

rule create_dataset_symlinks:
    output:
        os.path.join(out_dir, "input", "{dataset}.feather")
    run:
        os.symlink(dset_mapping[wildcards["dataset"]], output[0])
