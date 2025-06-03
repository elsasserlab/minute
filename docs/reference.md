# Reference guide

## Use of minute command-line interface

    minute <command> <options> <snakemake_specific_options>

**`command`** can be one of:

- **`init`**. Initializes a minute run creating the basic file structure and
template configurations.
- **`download`**. Downloads FASTQ files from SRA.
- **`run`**. Runs the workflow.

## General options

**`-h`**, **`--help`** Show help and exit.

**`--version`** Show version number and exit.


## minute commands

### init

Initialize a minute run. Parameters:

#### Positional

**`first`**. Path to the target directory where minute will run.

#### Named

**`--reads`** Path to the directory where the FASTQ file pairs are located.

**`--barcodes`** Path to the barcodes.tsv file. If not provided, `init` will
not create libraries.tsv, groups.tsv files. Requires `--reads` to be provided.

**`--input`** Prefix of the FASTQ file pair that corresponds to the Input sample.
It must match a file pair in the path specified by `--reads`. Required if 
`--barcodes` is provded.  

**`--config`** Path to a `minute.yaml` file to reuse. It will be copied into the
target directory.


### download

Download FASTQ files from the Sequence Read Archive (SRA). It will download
files with SRR codes specified in the fourth column of libraries.tsv file. 

Parameters: None.

### run

Run an already configured minute experiment directory. This is a wrapped 
Snakemake call, and its parameters work as Snakemake parameters. For instance,
for a minute dry run of the target `full`:

    minute run --dry-run full

See [Snakemake documentation](https://snakemake.readthedocs.io/) for more details.


#### Alternative target rules

By default, `minute run` will run the `final` rule, which generates demultiplexed
FASTQ file pairs, final BAM files and scaled bigWig files, plus unscaled (1x 
genome coverage RPGC) for comparison. Note that Input bigWig tracks are always
unscaled.

Alternative to this are:

**`full`**. Includes computationally costly extra QC: deepTools fingerprint.

**`quick`**. Runs the default pipeline but reduces the amount of bigWig files
produced. Only scaled bigWig files will be generated for treatments, and
unscaled bigWig files will be generated for controls.

**`no_bigwigs`**. Runs demultiplexing and alignment, but no bigWig generation. This
still produces full scaling information and MultiQC report, and it is very quick
to run.

**`mapq_bigwigs`**. Produces extra bigWigs `mapq.bw` filtered by mapping quality.
This requires a `mapping_quality_bigwig` parameter value higher than zero in the
`minute.yaml` configuration.

**`pooled_only`**. Produces only the bigWigs for the pooled replicates.

**`pooled_only_minimal`**. Runs the default pipeline, but it only generates
bigWig files for the pools present in groups.tsv, scaled for treatment libraries
and unscaled for the controls.

**`no_bigwigs`**. Runs the default pipeline but skips the bigWig generation
entirely. This is useful for a quick check on QC metrics, since the rest of the
steps are significantly faster than bigWig generation.

**`no_scaling`**. Skips entirely the scaling statistics and bigWigs. Produces
only unscaled bigWigs and a MultiQC report with only non-scaling QC metrics.
