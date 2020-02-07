#!/bin/bash
set -euo pipefail

rm -r testrun
mkdir testrun
cd testrun
cp ../minute-testdata/* .
mkdir fastq
mv *.fastq.gz fastq/
bowtie2-build ref.fa.gz ref
snakemake -j 1 -s ../Snakefile
