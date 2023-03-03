# minute Changelog

## v0.4.0

### Features

* PR #164: Bowtie2 indexes are now explicit on `minute.yaml` file as a
`bowtie2_index` field, and will be dynamically generated if this field is
missing or the path provided is not valid.

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
