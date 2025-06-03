# minute Changelog

## v0.11.0

### Features

* Introduced pooled_only_minimal rule that generates only scaled replicate pool
for treatments, and unscaled replicate pools for controls in the groups.tsv
file.
* Rule quick now produces the equivalent minimal set of bigWig (independent 
replicate files: scaled for treatment, unscaled for controls).

## v0.10.1

### Bug fixes

* Fix edge-case exception when trying to convert NA value to int after parsing
an empty insert size metrics file.

### Other

* Prioritize reporting rules over bigWig generation rules.

## v0.10.0

### Features

* New target rule `no_scaling` to skip the scaling process completely and
produce only the unscaled bigWig files and a final MultiQC report with
non-scaling related QC metrics.

## v0.9.2

### Minor
* More thorough check of user input, including: validating minute.yaml
has required fields, reference files exist, FASTQ files match libraries.tsv,
bowtie2 index is not required if aligner strobealign is used.
* More verbose command line messages when an error occurs.
* More verbose docstring for `minute run` command.

## v0.9.1

### Features
* Handle barcode dropouts. Zero values appear in the final report, and
this does not prevent final stats report to be produced.
* Carry barcode sequence to stats_summary.txt

### Other

* Minor: Styling changes made to the final Minute barcode barplots

## v0.9.0

### Features

* Bowtie2 run mode can be configured in `minute.yaml` to any of the bowtie2 
accepted modes. Default is still `--fast`.
* Optionally align reads with `strobealign` instead of `bowtie2`.
* Add read group RG to output BAM headers and read group ID to individual
alignments.

## v0.8.0

### Features

* New target rule `pooled_only` to skip individual replicates and generate 
bigwig tracks only for replicate pools.

## v0.7.0

### Features

* `minute init` now can accept a `--config` parameter to copy an already
existing `minute.yaml` into the run directory.

### Bug fixes

* Improved parameter check of `minute init`. Now a directory will not be created
if an error occurs.

## v0.6.0

### Features

* Barcode representation barplots, final number of mapped reads and percentage
of total.
* Stats summary table is now separated in replicates and pool for easier reading.

### Other
* Minor: expanded library name column so it is now wider.
* Automatic formatting of QC pass for insert size and duplication rates are 
removed.

## v0.5.0

### Features

* Added `mapping_quality` optional parameter to filter out reads under a certain
mapping quality threshold. By default, this value is set to zero (keep all).
* Added `mapping_quality_bigwig` specific to `mapq_bigwigs` new rule that
generates both final unfiltered bigWigs and `mapq.bw` filtered bigWigs according
to this parameter. By default, those filtered bigWigs are not generated.

## v0.4.1

### Bug fixes

* Handle missing/incorrect `--input` parameter values on `minute init`.

## v0.4.0

### Features

* Added `--barcodes` and `--input` parameters to `minute init` functionality to
streamline configuration of standard `minute` runs.

## v0.3.1

### Bug fixes

* Issue #159: Temporary removed dynamically renaming of dynamically generated
rules to avoid issues with recent versions of snakemake.

## v0.3.0

### Bug fixes

* Recover missing rows from Statistics Summary table on the MultiQC report when
a sample is mapped to more than one reference.
* Issue #159: Pin snakemake version to 7.16. 

### Features

* Improved MultiQC reporting by adding:
	- Picard `CollectInsertSizeMetrics` figure.
	- Custom library complexity estimate based on Picard, adapted to use the
	single-end numbers from Je.
	- Number of demultiplexed reads per sample.
	- Reorder tables on MultiQC report (Statistics Summary first).

### Other

* Minor numerical formatting changes on final stats summary tables.
