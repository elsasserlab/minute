# Minute pipeline

## Setup

- Install Conda and Bioconda. Use the
  [Bioconda instructions](https://bioconda.github.io/user/install.html) if you
  don’t have Conda and/or Bioconda, yet.
- Create a Conda environment with the necessary programs:

        conda env create -n minute -f environment.yml

- Activate the environment

        conda activate minute


## Running

1. Create an empty folder somewhere (called `myexperiment` in the following)
2. Create a subfolder `fastq` in the `myexperiment` folder
3. Create symbolic links within the `fastq/` folder that point to your input
   FASTQ files (those that contain the multiplexed libraries). Read 1 and
   read 2 files must have names ending in `_R1.fastq.gz` and `_R2.fastq.gz`,
   respectively.
4. Copy the `config.yaml` file into `myexperiment/` and edit it as required.
5. Create an `experiment.tsv` file describing your libraries (see below).
6. Run `snakemake -s path/to/the/Snakefile`.


## The experiment.tsv file

`experiment.tsv` is a text file in tab-separated value format describing the
sequenced libraries, one row per library. It needs to be placed in the
`myexperiment/` directory.

The columns are:

1. Sample name
2. Replicate id
3. Barcode
4. FASTQ base name

Use the provided `experiment.tsv` example file and edit it.

The FASTQ base name refers to files within the `fastq/` folder. The suffixes
`_R1.fastq.gz` and `_R2.fastq.gz` will be added automatically.


## The config.yaml file

The `config.yaml` is used to configure everything else. Create a copy of the
example `config.yaml` and edit it to suit your needs.
