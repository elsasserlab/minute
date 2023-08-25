#!/bin/bash
set -euo pipefail

version="0.10"

python -m pytest tests/

testdata_file=minute-testdata-${version}.tar.gz
if ! [[ -f ${testdata_file} ]]; then
  wget -nv https://export.uppmax.uu.se/snic2020-6-3/minute-testdata/${testdata_file}
fi

rm -rf testrun
minute init testrun
cd testrun
tar --strip-components=1 -xf ../${testdata_file}
cp ../testdata/*.{tsv,bed,sizes,yaml} .
mkdir fastq && mv -- *.fastq.gz fastq/
mkdir ref && mv ref*.fa.gz ref/
( cd ref && bowtie2-build -q ref1.fa.gz ref1 && bowtie2-build -q ref2.fa.gz ref2 )

minute run "$@"

# Run a standard test
cd ..
rm -rf data_std
rm -rf testrun_std

mkdir data_std && cd data_std
tar --strip-components=1 -xf ../${testdata_file}
mkdir fastq && mv -- testdata-*.fastq.gz fastq/

cd ..
minute init testrun_std --reads ./data_std/fastq --barcodes ./testdata/barcodes.tsv --input testdata-Inp
cd testrun_std
cp ../testdata/*.{bed,sizes,yaml} .
mkdir ref && mv ../data_std/ref*.fa.gz ref/
( cd ref && bowtie2-build -q ref1.fa.gz ref1 && bowtie2-build -q ref2.fa.gz ref2 )
minute run "$@"
