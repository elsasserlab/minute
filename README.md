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
5. Create `libraries.tsv` and `groups.tsv` files describing your libraries
   (see below).
6. Run `snakemake -p -s path/to/the/Snakefile`.

The `-p` makes `snakemake` print the commands that it executes and can be
omitted if you don’t want to see that.

Additionally, you may want to run `snakemake` with the option `-n` first, which
will only show which steps would be executed and not actually run them
(“dry run”).


## Configuration files

The configuration files `libraries.tsv`, `groups.tsv` and `config.yaml`
need to be placed in the `myexperiment/` directory. Use the provided
example files as templates and adjust as needed.

### The libraries.tsv file

`libraries.tsv` is a text file in tab-separated value format describing the
sequenced libraries, one row per library.

The columns are:

1. Sample name
2. Replicate id
3. Barcode
4. FASTQ base name

The FASTQ base name refers to files within the `fastq/` folder. The suffixes
`_R1.fastq.gz` and `_R2.fastq.gz` will be added automatically.


## The groups.tsv file

The `groups.tsv` file is a text file in tab-separated value format that
defines which of the libraries are the treatments and which are the
controls (or “inputs”).

The columns are are:

1. Treatment sample name
2. Replicate id
3. Control sample name


## The config.yaml file

The `config.yaml` is used to configure everything else. Please open the
file in an editor, read through the comments and edit as required.


## Running on SLURM clusters (Uppmax in Sweden)

Snakemake supports running on HPC environments. As such, it is possible to
run minute on SLURM clusters, including the Swedish UPPMAX clusters. Handling
of `config.yaml` and `libraries.tsv` files
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
- Make sure you call `conda activate minute` in the wrapper.


# Result folders and files

* `reports/multiqc_report.html`: MultiQC report
* `final/bam`: Final BAM files
* `final/bigwig`: Scaled and unscaled BigWig files
* `final/fastq`: Demultiplexed FASTQ files.
* `tmp/`: Intermediate files.
