# minute Changelog

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
