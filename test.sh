#!/bin/bash
set -euo pipefail

rm -rf testrun
mkdir testrun
cd testrun
cp ../minute-testdata/* .
mkdir fastq
mv *.fastq.gz fastq/
bowtie2-build ref.fa.gz ref
snakemake -p -j 1 -s ../Snakefile "$@"
