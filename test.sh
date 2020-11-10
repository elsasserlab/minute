#!/bin/bash
set -euo pipefail

version="0.8"

testdata_file=minute-testdata-${version}.tar.gz
if ! [[ -f ${testdata_file} ]]; then
  wget -nv https://export.uppmax.uu.se/snic2020-6-3/minute-testdata/${testdata_file}
fi

rm -rf testrun
mkdir testrun
cd testrun
tar --strip-components=1 -xf ../${testdata_file}
cp ../testdata/*.{tsv,bed,sizes,yaml} .
mkdir fastq && mv -- *.fastq.gz fastq/
mkdir ref && mv ref.* ref/
( cd ref && bowtie2-build -q ref.fa.gz ref )

snakemake -p -j 2 -s ../Snakefile "$@"
