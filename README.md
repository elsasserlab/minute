# Minute pipeline

## Setup

- Install Conda
- Create a Conda environment with the necessary programs:

        conda env create -n minute -f environment.yml

- Activate the environment

        conda activate minute

## Input

- `experiment.tsv` is a text file in tab-separated value format describing the
  sequenced libraries, one row per library. Columns:
    1. Sample name
    2. Replicated id
    3. Barcode
    4. FASTQ base name
- `config.yaml` is for all other paths and settings. Example:

        indexed_reference: "ref/mm9"
        blacklist_bed: "blacklist.bed"
