#!/bin/bash
set -euo pipefail

python -m pytest tests/

rm -rf testrun
minute init testrun
cd testrun
cp ../testdata/*.{tsv,bed,yaml} .
mkdir fastq && cp ../testdata/fastq/*.fastq.gz fastq/
mkdir ref && cp ../testdata/ref/ref*.fa.gz ref/
( cd ref && bowtie2-build -q ref1.fa.gz ref1 && bowtie2-build -q ref2.fa.gz ref2 )

minute run "$@"

# Run a standard test
cd ..
rm -rf data_std
rm -rf testrun_std
mkdir data_std && cd data_std
mkdir fastq && cp ../testdata/fastq/testdata-*.fastq.gz fastq/
cd ..
minute init testrun_std --reads data_std/fastq --barcodes testdata/barcodes.tsv --input testdata-Inp --config testdata/minute.yaml
cd testrun_std
cp ../testdata/*.bed .
mkdir ref && cp ../testdata/ref/ref*.fa.gz ref/
( cd ref && bowtie2-build -q ref1.fa.gz ref1 && bowtie2-build -q ref2.fa.gz ref2 )

minute run "$@"
