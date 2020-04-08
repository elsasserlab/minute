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
6. Run `snakemake -p -s path/to/the/Snakefile`.

The `-p` makes `snakemake` print the commands that it executes and can be
omitted if you don’t want to see that.

Additionally, you may want to run `snakemake` with the option `-n` first, which
will only show which steps would be executed and not actually run them
(“dry run”).


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

## Running on SLURM clusters (Uppmax in Sweden)

Snakemake supports running on HPC environments. As such, it is possible to
run minute on SLURM clusters, including the Swedish UPPMAX clusters. Handling
of `config.yaml` and `experiment.tsv` files
will work the same. You just need to have conda available
and an active minute environment that you can install as described in the
**setup** section.

**Note:** creating a Conda environment on a shared filesystem can be a
lengthy process.
It is not advised to run this on a login node. You can wrap your conda
environment installation on a `sbatch` file and submit it to the queue.

To specify cluster-specific parameters you can create a `cluster.yaml`
config file with your defaults:

        __default__:
            partition: "core"
            cpus: "{threads}"
            time: "0-06:00:00"
            project: "snicYYYY-NNN-N"
            jobname: "{rule}_{jobid}"
            
Note that you can use rule-dependent parameters such as `rule`, `jobid` and
`threads`. Then you call `snakemake`:

        snakemake -p -s path/to/the/Snakefile --jobs 20 --cluster-config path/to/cluster.yaml --cluster 'sbatch -A {cluster.project} -t {cluster.time} -c {cluster.cpus} -e logs_slurm/{cluster.jobname}.err -o logs_slurm/{cluster.jobname}.out -J {cluster.jobname}'

The `project` field is required, as SLURM will not queue your jobs if they are
not attached to a computing project. The `--jobs` parameter in the `snakemake`
command limits the maximum number of jobs to be queued, and it's also required.

In this example, separate log files for each job (stdout and stderr) are written
to a `logs_slurm` folder (otherwise you get a bunch of `slurm-<jobid>.out` files
in the working directory). If you want this behavior, you need to create that
directory before running `snakemake`.

You can also wrap your minute pipeline call in a `sbatch` file itself, so
the scheduler does not run on a login node:

- Create the `logs_slurm` directory in the `myexperiment` directory.
- Assign a `-t` (time) parameter long enough, as this will be mostly waiting
  for other jobs to run. Consider that this includes not only running time but
  also queuing time.
- Ask only for 1 core.
- Make sure you call `conda init && conda activate minute` in the wrapper.


# Result folders

* `fastqc/`: FastQC report for input reads
* `noumi/`: Input reads without UMIs
* `demultiplexed/`: Demultiplexed “noumi” reads
* `mapped/`: Demultiplexed reads mapped to reference
* `mapped_se/`: Mapped reads converted to single end
* `dupmarked/`: Duplicate-marked (single-end) reads
* `restricted/`: Mapped reads restricted to regions of interest. Currently this
  folder is the one that contains the final BAMs.
