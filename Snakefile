# TODO
# - switch to interleaved files?
from itertools import groupby
from utils import read_experiment_description


configfile: "config.yaml"

libraries = list(read_experiment_description())

rule all:
    input:
        expand([
            "igv/{library.name}.bw",
            "igv/{library.name}.tdf",
            "restricted/{library.name}.bam",
            "results/{library.name}.insertsizes.pdf",
            "results/{library.name}.insertsizes.txt",
        ], library=libraries)


# TODO this needs to be replaced with a proper tool
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
    """File with list of barcodes needed for demultiplexing"""
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
        bam="mapped/{library}.bam"
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


# TODO have a look at UMI-tools also
rule mark_duplicates:
    """UMI-aware duplicate removal with je suite"""
    output:
        bam="dupmarked/{library}.bam",
        metrics="dupmarked/{library}.metrics"
    input:
        bam="mapped/{library}.bam"
    shell:
        """
        je markdupes MISMATCHES=1 REMOVE_DUPLICATES=TRUE SLOTS=-1 I={input.bam} O={output.bam} M={output.metrics}
        """


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
    shell:
        "bamCoverage"
        " -p max"
        " --normalizeUsing RPGC"
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
