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

rule all:
    input:
        expand(os.path.join(out_dir, "results/networks/{dataset}.graphml"),
                dataset=dset_names)

rule summarize_results:
    input:
        os.path.join(out_dir, "results/merged/{dataset}.feather")
    output:
        summary=os.path.join(out_dir, "results/summary/{dataset}_summary.tsv"),
    script: "src/summarize_results.R"

rule create_networks:
    input:
        os.path.join(out_dir, "results/{dataset}/gene_sets.feather")
    output:
        os.path.join(out_dir, "results/networks/{dataset}.graphml")
    script: "src/create_networks.R"

rule run_fgsea:
    input: 
        dataset=os.path.join(out_dir, "input", "{dataset}.feather"),
    output:
        os.path.join(out_dir, "results/{dataset}/gene_sets.feather")
    script: "src/run_fgsea.R"

rule create_dataset_symlinks:
    output:
        os.path.join(out_dir, "input", "{dataset}.feather")
    run:
        os.symlink(dset_mapping[wildcards["dataset"]], output[0])
