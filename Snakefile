import se_bam
# TODO
# - switch to interleaved files?
from itertools import groupby
from utils import (
    read_libraries,
    read_controls,
    flagstat_mapped_reads,
    compute_scaling,
    parse_duplication_metrics,
    parse_insert_size_metrics,
)


configfile: "config.yaml"

libraries = list(read_libraries())
normalization_pairs = list(read_controls(libraries))  # or: normalization_groups

# Map a FASTQ prefix to its list of libraries
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
            "stats/{library.name}.txt",
        ], library=libraries),
        expand("fastqc/{fastq}_R{read}_fastqc.html",
            fastq=fastq_map.keys(), read=(1, 2)),
        expand("scaled/{library.name}.scaled.bw",
            library=[np.treatment for np in normalization_pairs]),


rule clean:
    shell:
        "rm -rf"
        " barcodes"
        " noumi"
        " demultiplexed"
        " mapped"
        " mapped_se"
        " dupmarked"
        " dedup"
        " results"
        " restricted"
        " igv"
        " fastqc"
        " scaled"
        " stats"
        " factors"
        " log"


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
            fastqbase=fastq_base,
        log:
            "log/demultiplexed/{fastqbase}.log".format(fastqbase=fastq_base)
        shell:
            "cutadapt"
            " -e 0.15"  # TODO determine from barcode length
            " -g file:{input.barcodes_fasta}"
            " -o {params.r1}"
            " -p {params.r2}"
            " --untrimmed-output demultiplexed/{params.fastqbase}-unknown_R1.fastq.gz"
            " --untrimmed-paired-output demultiplexed/{params.fastqbase}-unknown_R2.fastq.gz"
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
        r2="demultiplexed/{library}_R2.fastq.gz",
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


rule convert_to_single_end:
    """Convert sam files to single-end for marking duplicates"""
    output:
        bam="mapped_se/{library}.bam"
    input:
        bam="mapped/{library}.bam"
    run:
        se_bam.convert_paired_end_to_single_end_bam(input.bam, output.bam)


# TODO have a look at UMI-tools also
rule mark_duplicates:
    """UMI-aware duplicate marking with je suite"""
    output:
        bam="dupmarked/{library}.bam",
        metrics="dupmarked/{library}.metrics"
    input:
        bam="mapped_se/{library}.bam"
    shell:
        "LC_ALL=C je"
        " markdupes"
        " MISMATCHES=1"
        " REMOVE_DUPLICATES=FALSE"
        " SLOTS=-1"
        " SPLIT_CHAR=_"
        " I={input.bam}"
        " O={output.bam}"
        " M={output.metrics}"


rule deduplicate_pe_file:
    """Select duplicate-flagged alignments and filter in the PE file"""
    output:
        bam="dedup/{library}.bam"
    input:
        target_bam="mapped/{library}.bam",
        proxy_bam="dupmarked/{library}.bam"
    run:
        se_bam.mark_duplicates_by_proxy_bam(
            input.target_bam,
            input.proxy_bam,
            output.bam)


rule remove_exclude_regions:
    output:
        bam="restricted/{library}.bam"
    input:
        bam="dedup/{library}.bam",
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
        txt="restricted/{library}.insertsizes.txt",
        pdf="restricted/{library}.insertsizes.pdf",
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
rule bigwig:
    output:
        bw="igv/{library}.bw"
    input:
        bam="restricted/{library}.bam",
        bai="restricted/{library}.bai"
    threads: 1
    shell:
        "bamCoverage"
        " -p {threads}"
        " --normalizeUsing RPGC"
        " --effectiveGenomeSize {config[genome_size]}"
        " -b {input.bam}"
        " -o {output.bw}"


rule compute_scaling_factors:
    input:
        treatments=["restricted/{library.name}.flagstat.txt".format(library=np.treatment) for np in normalization_pairs],
        controls=["restricted/{library.name}.flagstat.txt".format(library=np.control) for np in normalization_pairs],
    output:
        factors=["factors/{library.name}.factor.txt".format(library=np.treatment) for np in normalization_pairs],
        info="scalinginfo.txt"
    run:
        with open(output.info, "w") as outf:
            factors = compute_scaling(
                normalization_pairs,
                input.treatments,
                input.controls,
                outf,
                genome_size=config["genome_size"],
                fragment_size=config["fragment_size"],
            )
            for factor, factor_path in zip(factors, output.factors):
                with open(factor_path, "w") as f:
                    print(factor, file=f)


rule scaled_bigwig:
    output:
        bw="scaled/{library}.scaled.bw"
    input:
        factor="factors/{library}.factor.txt",
        bam="restricted/{library}.bam",
        bai="restricted/{library}.bai",
    threads: 4
    shell:
        # TODO also run this
        # - with "--binSize 50 --smoothLength 150"
        # - with "--binSize 500 --smoothLength 5000"
        "bamCoverage"
        " -p {threads}"
        " --binSize 1"
        " --extendReads"
        " --scaleFactor $(< {input.factor})"
        " --bam {input.bam}"
        " -o {output.bw}"


rule stats:
    output:
        txt="stats/{library}.txt"
    input:
        mapped="mapped/{library}.flagstat.txt",
        dupmarked="dupmarked/{library}.flagstat.txt",
        restricted="restricted/{library}.flagstat.txt",
        metrics="dupmarked/{library}.metrics",
        insertsizes="restricted/{library}.insertsizes.txt",
    run:
        row = []
        for flagstat, name in [
            (input.mapped, "mapped"),
            (input.dupmarked, "dupmarked"),
            (input.restricted, "restricted"),
        ]:
            mapped_reads = flagstat_mapped_reads(flagstat)
            row.append(mapped_reads)
        row.append(parse_duplication_metrics(input.metrics)["estimated_library_size"])
        row.append(parse_insert_size_metrics(input.insertsizes)["median_insert_size"])
        with open(output.txt, "w") as f:
            print("mapped", "dupmarked_mapped", "restricted_mapped", "library_size", "insert_size", sep="\t", file=f)
            print(*row, sep="\t", file=f)


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
