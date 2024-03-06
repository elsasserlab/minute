# User guide

Broadly, the input of a **minute** run consists of a set of multiplexed paired-end
FASTQ files which are demultiplexed and aligned to a reference genome (with
quality-control metrics and intermediate steps such as deduplication, optional
filtering by mapping quality and/or user-defined excluded loci). The final
numbers of reads aligned for both the input and each ChIP are used to normalize
and produce MINUTE-scaled bigWig files. 


## Standard run

This is the most comon use-case for **minute** users, and the recommended 
starting point.

### Configuration

In order to process your data, **`minute`** needs the barcode sequence corresponding
to each library in the pool, and which of them is used as reference for scaling.
You can simply specify these in a tab-separated `barcodes.tsv` file:

    # Barcode information file
    # Columns: barcode name, replicate id, barcode, reference
    reference_condition  1       CCGGCGTT        mm39
    reference_condition  2       TTATTTCT        mm39
    treatment_condition  1       GATGTTTG        mm39
    treatment_condition  2       TTGTGAAA        mm39


This file has the following columns:

1. **Condition name**. Any name you want to attach to the barcode. It should not 
contain space characters.
2. **Replicate id**. Different replicates for the same condition name.
3. **Barcode sequence**. Nucleotide sequence that defines the pair condition, replicate.
4. **Genome reference**. Reference identifier to map the reads to. This identifier should match an 
identifier in the `minute.yaml` [file](guide.md#the-minuteyaml-file). 

#### Configuration steps

1. Create a barcodes `barcodes.tsv` file describing the barcodes for each 
condition of your experiment as described.
2. Activate `minute` conda environment (see the [installation guide](install.md)).

         conda activate minute

3. Run:

         minute init myexperiment --reads path/to/fastq --barcodes path/to/barcodes.tsv --input inp_prefix
Note that `inp_prefix` **must match one of the FASTQ read pairs** in your `fastq` directory: `inp_prefix_R1.fastq.gz`, `inp_prefix_R2.fastq.gz`.

      This will create a `myexperiment` directory that contains a `fastq` subdirectory
with symlinks to each FASTQ file in `path/to/fastq`, a template `minute.yaml`
configuration file and `libraries.tsv`, `groups.tsv` files.

      Optionally, if you have ran `minute` before in your system and you already
have some `minute.yaml` configuration you want to reuse, you can specify it via
the `--config` parameter. This will make a copy of the provided file in the 
directory where `minute` will run.

         minute init myexperiment \
            --reads path/to/fastq \
            --barcodes path/to/barcodes.tsv \
            --input inp_prefix \
            --config path/to/minute.yaml

4. Edit [`minute.yaml`](guide.md#the-minuteyaml-file) if required. References present in `barcodes.tsv` must be 
specified in it. For example, if you specify `mm39` as reference genome, there 
must be a `mm39` entry in `minute.yaml` with paths to matching bowtie2 indexes
and FASTA reference. 


!!! Note

    **minute** runs on paired-end FASTQ files. Read 1 and read 2 files must have names ending in `_R1.fastq.gz` and `_R2.fastq.gz`,
    respectively.

This configuration is enough for a standard `minute` run.


### Running minute

Once you have set up a `minute` run, your experiment directory should look
similar to this (with your own FASTQ file pairs):

    myexperiment
      ├── fastq
      │   ├── H3K4me3-ChIP_R1.fastq.gz
      │   ├── H3K4me3-ChIP_R2.fastq.gz
      │   ├── INPUT_R1.fastq.gz
      │   └── INPUT_R2.fastq.gz
      ├── groups.tsv
      ├── libraries.tsv
      └── minute.yaml


At this point, you just have to move to your experiment directory and run minute:

    cd myexperiment && minute run


!!! Note

    You probably want to run `minute run` with the option `-n` first, which
    will do a “dry run”, that is, it only shows which steps would be executed and
    does not actually run them.


### Result folders and files

* `reports/multiqc_report.html`: MultiQC report
* `final/bam`: Final BAM files
* `final/bigwig`: Scaled and unscaled BigWig files.
* `final/fastq`: Demultiplexed FASTQ files.
* `tmp/`: Intermediate files.


## Advanced options

It is possible to further customize a `minute` run. You would only need to do
this if you have any of these following specific use cases:

1. Scaling groups do not correspond 1:1 to each FASTQ R1/R2 pair.
2. Some FASTQ R1/R2 pairs have different barcode/conditions.
3. The reference is not always the same on each scaling group, and/or you do not
want to use the pooled sample final read count as scaling factor.
4. You would like to align some samples to a different reference.

Note that in any of these cases it is probably easier to further edit the
configuration obtained from a `minute init` run,
so it is still recommended to first configure a [standard run](guide.md#standard-run).


### Custom configuration

Alternatively, you can edit or manually create the `libraries.tsv` and
`groups.tsv` as described. 

1. Create an empty folder somewhere (called `myexperiment` in the following)
2. Create a subfolder `fastq` in the `myexperiment` folder
3. Create symbolic links within the `fastq/` folder that point to your input
   FASTQ files (those that contain the multiplexed libraries). 
4. Copy the `testdata/minute.yaml` file into `myexperiment/` and edit it as required.
5. Create `libraries.tsv` and `groups.tsv` files describing your libraries
   (see below).

Once the custom configuration is done, you can proceed to [run it](guide.md#running-minute).
You can see more details about the configuration files below.


### Configuration files

The configuration files `libraries.tsv`, `groups.tsv` and `minute.yaml`
need to be placed in the `myexperiment/` directory. Use the provided
example files under `testdata` as templates and adjust as needed.

#### The libraries.tsv file

This is a text file in tab-separated value format describing the
sequenced libraries, one row per library.


    # Columns: sample name, replicate id, barcode, FASTQ base name
    #
    # If the barcode is ".", the FASTQ file will be used directly.
    # If the barcode is set to a nucleotide sequence, the FASTQ file
    # will be demultiplexed.
    #
    H3K4me3_reference_condition   1  CCGGCGTT H3K4me3-ChIP
    H3K4me3_reference_condition   2  TTATTTCT H3K4me3-ChIP
    H3K4me3_treatment_condition   1  GATGTTTG H3K4me3-ChIP
    H3K4me3_treatment_condition   2  TTGTGAAA H3K4me3-ChIP
    Input_reference_condition   1  CCGGCGTT INPUT
    Input_reference_condition   2  TTATTTCT INPUT
    Input_treatment_condition   1  GATGTTTG INPUT
    Input_treatment_condition   2  TTGTGAAA INPUT


The columns are:

1. **Sample name**. Any name you want to give to the sample. 
2. **Replicate number**.
3. **Barcode sequence** if applicable, a dot (`.`) otherwise. If a `.` is specified, the 
demultiplexing step of this library will be skipped. This is a feature that
allows to reprocess already demultiplexed MINUTE-ChIP datasets.
4. **FASTQ base name**. This must match the prefix of the FASTQ file pairs that
are going to be processed (without `_R1.fastq.gz`, `_R2.fastq.gz`).

Sample names in the first column must include some
information about the matched FASTQ file pair, since no two pairs sample_name
plus replicate can be repeated in the samplesheet. As this is a multiplexed
experiment, it is usually the case that condition+replicate+barcode are the
same on each ChIP and Input. If that is the case,
configuring your run using the [standard run](guide.md#standard-run) will
simplify the process.

The FASTQ base name refers to files within the `fastq/` folder. The suffixes
`_R1.fastq.gz` and `_R2.fastq.gz` will be added automatically.

If the barcode is ".", the FASTQ file will be used directly.
If the barcode is set to a nucleotide sequence, the FASTQ file pair
will be demultiplexed.

In the example above, the expected FASTQ files according to the fourth column
are `H3K4me3-ChIP_R1.fastq.gz`, `H3K4me3-ChIP_R2.fastq.gz`, `INPUT_R1.fastq.gz`, `INPUT_R2.fastq.gz`

If you use an SRA accession for the FASTQ base name, a separate command can be
used to automatically download the listed samples from the SRA, see [downloading data from SRA](guide.md#downloading-data-from-the-sequence-read-archive-sra).


#### The groups.tsv file

This is a text file in tab-separated value format that
defines which of the libraries are the treatments, which are the
controls (or “inputs”) and to which reference they need to be scaled to.
It is possible to specify different scaling groups, each one with its
own reference library. For each scaling group, the first library will
be taken as reference. Note that a given library cannot be specified in
multiple scaling groups at the same time. 

    # Columns: treatment name, replicate id, control name, scaling group, ref genome
    # 
    H3K4me3_reference_condition   pooled  Input_reference_condition       H3K4me3  mm39
    H3K4me3_treatment_condition   pooled  Input_treatment_condition       H3K4me3  mm39


The columns are:

1. **Treatment sample name**.
2. **Replicate id** - This column can also contain the reserved word `pooled`
to indicate that the replicates are to be pooled before scaling. 
3. **Control sample name**.
4. **Scaling group**.
5. **Reference genome**. Reads will be mapped to this genome. This identifier
must match a `reference` in the `minute.yaml` file. 


#### The minute.yaml file

The `minute.yaml` is used to configure everything else. If you initialized your
experiment as a Standard run using `minute init` you should have a template 
`minute.yaml` file on your `myexperiment` directory. Please open the file in
an editor, read through the comments and edit as required.

      # Configuration settings for the minute pipeline
      # Paths (unless absolute) are relative to the directory in which snakemake is run.
      references:
        mini:  # Arbitrary name for this reference. This is also used in output file names.

          # Path to a reference FASTA file (may be gzip-compressed).
          # A matching Bowtie2 index must exist in the same location.
          fasta: "ref/ref1.fa.gz"

          # Path to a BED file with regions to exclude.
          # If this is empty, no regions are excluded.
          exclude: "exclude1.bed"

        small:
          fasta: "ref/ref2.fa.gz"
          exclude:


      # Length of the 5' UMI
      umi_length: 6

      # Fragment length (insert size)
      fragment_size: 150

      # Allow this many errors in the barcode when demultiplexing
      max_barcode_errors: 1

      # Filter out reads under certain mapping quality (0: no filtering)
      mapping_quality: 0

      # If filtered_bigwigs rule is run, additional bigWig files are produced,
      # where reads with MQ lower than this value are filtered out. This parameter
      # is independent of mapping_quality parameter above.
      mapping_quality_bigwig: 20


## Downloading data from the Sequence Read Archive (SRA)

To use data from SRA that you have not already downloaded, create the
`libraries.tsv` file listing your libraries as described above.
For all samples to be downloaded from the SRA, write a dot (".") in
the *barcode* (third) column and put the SRA run accession
(starting with SRR, ERR, or DRR) into the *FASTQ base name* (fourth) column.

The `minute.yaml` and `groups.tsv` files are not needed for this step.

Then run:

    minute download

This will download the files into the `fastq/` directory.

If your `libraries.tsv` refers to non-SRA entries and you have not
placed those files into the `fastq/` directory, you will get an error
message for them, but the other files will still be downloaded.

You can now move, copy or link the downloaded files into a separate folder
somewhere for later use. The next time they are needed, you can avoid the
download and instead use the existing files.

If you download the data yourself, note that Minute expects file names
ending in `_R1.fastq.gz` and `_R2.fastq.gz`. Because `fastq-dump`
creates files named `..._1.fastq.gz` and `..._2.fastq.gz`, they need to be
renamed.


## Running on SLURM clusters

Snakemake supports running on HPC environments. As such, it is possible to
run minute on [SLURM clusters](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#executing-on-slurm-clusters). Handling of `minute.yaml` and `libraries.tsv` files
will work the same. You just need to have conda available
and an active minute environment that you can install as described in the
**setup** section.

To specify cluster-specific parameters you can create a `cluster.yaml`
config file with your defaults:

        __default__:
            partition: "core"
            cpus: "{threads}"
            time: "0-06:00:00"
            project: "project-code"
            jobname: "{rule}_{jobid}"
            
Note that you can use rule-dependent parameters such as `rule`, `jobid` and
`threads`. Then you call `minute run`:

        minute run --jobs 20 --cluster-config path/to/cluster.yaml --cluster 'sbatch -A {cluster.project} -t {cluster.time} -c {cluster.cpus} -e logs_slurm/{cluster.jobname}.err -o logs_slurm/{cluster.jobname}.out -J {cluster.jobname}'

The `project` field is required, as SLURM will not queue your jobs if they are
not attached to a computing project. The `--jobs` parameter in the `snakemake`
command limits the maximum number of jobs to be queued, and it's also required.

In this example, separate log files for each job (stdout and stderr) are written
to a `logs_slurm` folder (otherwise you get a bunch of `slurm-<jobid>.out` files
in the working directory). If you want this behavior, you need to create that
directory before running `snakemake`.

You can also wrap your minute pipeline call in a `sbatch` file itself, so
the scheduler does not run on a login node:

1. Create the `logs_slurm` directory in the `myexperiment` directory.
2. Assign a `-t` (time) parameter long enough, as this will be mostly waiting
  for other jobs to run. Consider that this includes not only running time but
  also queuing time.
3. Ask only for 1 core.
4. Make sure you call `conda activate minute` in the wrapper.