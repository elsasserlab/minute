# kate: syntax Python;

# TODO
# - switch to interleaved files?

from typing import NamedTuple
from itertools import groupby


configfile: "config.yaml"



class Library(NamedTuple):
    sample: str
    replicate: str
    barcode: str
    fastqbase: str

    @property
    def name(self):
        return f"{self.sample}_replicate{self.replicate}"


def read_experiment_description():
    with open("experiment.tsv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            fields[1] = fields[1].lstrip("R")
            yield Library(*fields)


libraries = list(read_experiment_description())

demuxtargets = expand("demultiplexed/{library.name}_R{read}.fastq.gz", library=libraries, read=(1, 2))

bams = expand("bam/{library.name}.bam", library=libraries)

rule all:
    input:
        bams


rule move_umi_to_header:
    output:
        fastq="noumi/{fastqbase}.interleaved.fastq.gz"
    input:
        r1="fastq/{fastqbase}_R1.fastq.gz",
        r2="fastq/{fastqbase}_R2.fastq.gz",
    shell:
        "paste "
        " <(zcat {input.r1} | paste - - - -)"
        " <(zcat {input.r2} | paste - - - -)"
        " "
        "| awk -F $'\\t' -vOFS=$'\\n' "
        "'{{ umi=substr($2, 1, 6); "
        "sub(/\\s/, \":\" umi \" \", $1); "
        "sub(/\\s/, \":\" umi \" \", $5); "
        "$2 = substr($2, 7); "
        "$4 = substr($4, 7); }};1'"
        " | pigz > {output.fastq}"


rule barcodes:
    output:
        barcodes_fasta="barcodes/{fastqbase}.fasta"
    run:
        with open(output.barcodes_fasta, "w") as f:
            for library in libraries:
                if library.fastqbase != wildcards.fastqbase:
                    continue
                f.write(f">{library.name}\n^{library.barcode}\n")





for key, items in groupby(sorted(libraries, key=lambda lib: lib.fastqbase), key=lambda lib: lib.fastqbase):
    items = list(items)

    rule:
        output:
            expand("demultiplexed/{library.name}_R{read}.fastq.gz", library=items, read=(1, 2))
        input:
            fastq="noumi/{fastqbase}.interleaved.fastq.gz".format(fastqbase=key),
            barcodes_fasta="barcodes/{fastqbase}.fasta".format(fastqbase=key),
        params:
            r1=lambda wildcards: "demultiplexed/{name}_R1.fastq.gz",
            r2=lambda wildcards: "demultiplexed/{name}_R2.fastq.gz",
        log:
            "log/demultiplexed/{fastqbase}.log".format(fastqbase=key)
        shell:
            "cutadapt"
            " --interleaved"
            " -g file:{input.barcodes_fasta}"
            " -o {params.r1}"
            " -p {params.r2}"
            " {input.fastq}"
            " > {log}"


rule bowtie2:
    threads:
        20
    output:
        bam="bam/{library}.bam"
    input:
        r1="demultiplexed/{library}_R1.fastq.gz",
        r2="demultiplexed/{library}_R1.fastq.gz",
    log:
        "log/bowtie2-{library}.log"
    # TODO
    # - --sensitive (instead of --fast) woul be default
    # - write uncompressed BAM?
    # - filter unmapped reads directly? (samtools view -F 4 or bowtie2 --no-unal)
    # - add RG header
    shell:
        "bowtie2"
        " -p {threads}"
        " -x {config[indexed_reference]}"
        " -1 {input.r1}"
        " -2 {input.r2}"
        " --fast"
        " 2> {log}"
        " "
        "| samtools sort -o {output.bam} -"
