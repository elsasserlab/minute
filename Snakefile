# TODO
# - switch to interleaved files?
from itertools import groupby
from utils import read_experiment_description


configfile: "config.yaml"

libraries = list(read_experiment_description())

# Map a FASTQ prefix to list of libraries
fastq_map = {
    fastq_base: list(libs)
    for fastq_base, libs in
    groupby(sorted(libraries, key=lambda lib: lib.fastqbase), key=lambda lib: lib.fastqbase)
}

rule all:
    input:
        expand([
            "igv/{library.name}.bw",
            "igv/{library.name}.tdf",
            "restricted/{library.name}.bam",
            "restricted/{library.name}.flagstat.txt",
            "restricted/{library.name}.idxstats.txt",
            "results/{library.name}.insertsizes.pdf",
            "results/{library.name}.insertsizes.txt",
        ], library=libraries),
        expand("fastqc/{fastq}_R{read}_fastqc.html",
               fastq=fastq_map.keys(), read=(1, 2))


rule clean:
    shell:
        "rm -rf"
        " barcodes"
        " noumi"
        " demultiplexed"
        " mapped"
        " dupmarked"
        " results"
        " restricted"
        " igv"
        " fastqc"


rule fastqc_input:
    output:
        "fastqc/{name}_fastqc.html"
    input:
        fastq="fastq/{name}.fastq.gz"
    shell:
        "fastqc -o fastqc {input.fastq}"
        " && rm fastqc/{wildcards.name}_fastqc.zip"


rule move_umi_to_header:
    output:
        r1="noumi/{name}.1.fastq.gz",
        r2="noumi/{name}.2.fastq.gz",
    input:
        r1="fastq/{name}_R1.fastq.gz",
        r2="fastq/{name}_R2.fastq.gz",
    params:
        umistring="N" * config['umi_length']
    shell:
        "umi_tools "
        " extract"
        " --extract-method=string"
        " -p {params.umistring}"
        " -I {input.r1}"
        " --read2-in={input.r2}"
        " -S {output.r1}"
        " --read2-out {output.r2}"


rule barcodes:
    """File with list of barcodes needed for demultiplexing"""
    output:
        barcodes_fasta="barcodes/{fastqbase}.fasta"
    run:
        with open(output.barcodes_fasta, "w") as f:
            for library in libraries:
                if library.fastqbase != wildcards.fastqbase:
                    continue
                f.write(f">{library.name}\n^{library.barcode}\n")


for fastq_base, libs in fastq_map.items():

    rule:
        output:
            expand("demultiplexed/{library.name}_R{read}.fastq.gz",
                library=libs, read=(1, 2))
        input:
            r1="noumi/{fastqbase}.1.fastq.gz".format(fastqbase=fastq_base),
            r2="noumi/{fastqbase}.2.fastq.gz".format(fastqbase=fastq_base),
            barcodes_fasta="barcodes/{fastqbase}.fasta".format(fastqbase=fastq_base),
        params:
            r1=lambda wildcards: "demultiplexed/{name}_R1.fastq.gz",
            r2=lambda wildcards: "demultiplexed/{name}_R2.fastq.gz",
        log:
            "log/demultiplexed/{fastqbase}.log".format(fastqbase=fastq_base)
        shell:
            "cutadapt"
            " -g file:{input.barcodes_fasta}"
            " -o {params.r1}"
            " -p {params.r2}"
            " {input.r1}"
            " {input.r2}"
            " > {log}"


rule bowtie2:
    threads:
        20
    output:
        bam="mapped/{library}.bam"
    input:
        r1="demultiplexed/{library}_R1.fastq.gz",
        r2="demultiplexed/{library}_R1.fastq.gz",
    log:
        "log/bowtie2-{library}.log"
    # TODO
    # - --sensitive (instead of --fast) would be default
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


# TODO have a look at UMI-tools also
rule mark_duplicates:
    """UMI-aware duplicate removal with je suite"""
    output:
        bam="dupmarked/{library}.bam",
        metrics="dupmarked/{library}.metrics"
    input:
        bam="mapped/{library}.bam"
    shell:
        "je"
        " markdupes"
        " MISMATCHES=1"
        " REMOVE_DUPLICATES=TRUE"
        " SLOTS=-1"
        " SPLIT_CHAR=_"
        " I={input.bam}"
        " O={output.bam}"
        " M={output.metrics}"


rule remove_exclude_regions:
    output:
        bam="restricted/{library}.bam"
    input:
        bam="dupmarked/{library}.bam",
        bed=config["blacklist_bed"]
    shell:
        "bedtools"
        " intersect"
        " -v"
        " -abam {input.bam}"
        " -b {input.bed}"
        " > {output.bam}"


rule insert_size_metrics:
    output:
        txt="results/{library}.insertsizes.txt",
        pdf="results/{library}.insertsizes.pdf",
    input:
        bam="restricted/{library}.bam"
    shell:
        "picard"
        " CollectInsertSizeMetrics"
        " I={input.bam}"
        " O={output.txt}"
        " HISTOGRAM_FILE={output.pdf}"
        " MINIMUM_PCT=0.5"
        " STOP_AFTER=10000000"


rule igvtools_count:
    output:
        tdf="igv/{library}.tdf"
    input:
        bam="restricted/{library}.bam",
        chrom_sizes=config["chrom_sizes"]
    shell:
        "igvtools"
        " count"
        " --extFactor 60"
        " {input.bam}"
        " {output.tdf}"
        " {input.chrom_sizes}"


# TODO can genome_size be computed automatically?
# TODO this is slow on tiny test datasets
rule bigwig:
    output:
        bw="igv/{library}.bw"
    input:
        bam="restricted/{library}.bam",
        bai="restricted/{library}.bai"
    threads: 1
    shell:
        "bamCoverage"
        " -p {threads}" # TODO " --normalizeUsing RPGC"
        " --effectiveGenomeSize {config[genome_size]}"
        " -b {input.bam}"
        " -o {output.bw}"


rule samtools_index:
    output:
        "{name}.bai"
    input:
        "{name}.bam"
    shell:
        "samtools index {input} {output}"


rule samtools_idxstats:
    output:
        txt="{name}.idxstats.txt"
    input:
        bam="{name}.bam",
        bai="{name}.bai",
    shell:
        "samtools idxstats {input.bam} > {output.txt}"


rule samtools_flagstat:
    output:
        txt="{name}.flagstat.txt"
    input:
        bam="{name}.bam"
    shell:
        "samtools flagstat {input.bam} > {output.txt}"
