import os
import sys

import se_bam
from utils import (
    read_libraries,
    read_scaling_groups,
    flagstat_mapped_reads,
    compute_scaling,
    parse_duplication_metrics,
    parse_insert_size_metrics,
    parse_stats_fields,
    read_int_from_file,
    compute_genome_size,
    detect_bowtie_index_name,
    get_replicates,
    group_libraries_by_sample,
    format_metadata_overview,
    is_snakemake_calling_itself,
    map_fastq_prefix_to_list_of_libraries,
    get_normalization_pairs,
)


configfile: "config.yaml"


localrules: 
    final,
    multiqc,
    clean,
    barcodes,
    remove_exclude_regions,
    summarize_scaling_factors,
    extract_fragment_size,
    stats,
    stats_summary,
    convert_to_single_end,
    mark_duplicates,
    mark_pe_duplicates,
    samtools_index,
    samtools_flagstat,
    samtools_flagstat_final,
    insert_size_metrics,


if "bowtie_index_name" not in config:
    try:
        config["bowtie_index_name"] = detect_bowtie_index_name(config["reference_fasta"])
    except FileNotFoundError as e:
        sys.exit(str(e))

libraries = list(read_libraries())
pools = list(group_libraries_by_sample(libraries))
scaling_groups = list(read_scaling_groups(libraries))

if not is_snakemake_calling_itself():
    print(format_metadata_overview(libraries, pools, scaling_groups), file=sys.stderr)

fastq_map = map_fastq_prefix_to_list_of_libraries(libraries)


rule final:
    input:
        "reports/multiqc_report.html",
        expand("final/bigwig/{library.name}.unscaled.bw", library=libraries),
        expand("final/bigwig/{pool.name}.unscaled.bw", pool=pools),
        expand("final/bigwig/{library.name}.scaled.bw",
            library=[np.treatment for np in get_normalization_pairs(scaling_groups)]),


rule multiqc:
    output: "reports/multiqc_report.html"
    input:
        expand("reports/fastqc/{fastq}_R{read}_fastqc/fastqc_data.txt",
            fastq=fastq_map.keys(), read=(1, 2)),
        expand("log/2-noadapters/{fastq}.trimmed.log", fastq=fastq_map.keys()),
        expand("log/4-mapped/{library.name}.log", library=libraries),
        expand("stats/6-dupmarked/{library.name}.metrics", library=libraries),
        "reports/scalinginfo.txt",
        "reports/stats_summary.txt",
        multiqc_config=os.path.join(os.path.dirname(workflow.snakefile), "multiqc_config.yaml")
    shell:
        "multiqc -o reports/ -c {input.multiqc_config} {input}"


rule clean:
    shell:
        "rm -rf"
        " tmp"
        " final"
        " stats"
        " reports"
        " log"


rule fastqc_input:
    output:
        html="reports/fastqc/{name}_fastqc.html",
        zip=temp("reports/fastqc/{name}_fastqc.zip"),
        data="reports/fastqc/{name}_fastqc/fastqc_data.txt",
    input:
        fastq="fastq/{name}.fastq.gz"
    shell:
        "fastqc --extract -o reports/fastqc {input.fastq}"


rule move_umi_to_header:
    output:
        r1=temp("tmp/1-noumi/{name}.1.fastq.gz"),
        r2=temp("tmp/1-noumi/{name}.2.fastq.gz"),
    input:
        r1="fastq/{name}_R1.fastq.gz",
        r2="fastq/{name}_R2.fastq.gz",
    params:
        umistring="N" * config['umi_length']
    log:
        "log/1-noumi/{name}.log"
    group: "clean_reads"
    shell:
        "umi_tools"
        " extract"
        " --extract-method=string"
        " -p {params.umistring}"
        " -I {input.r1}"
        " --read2-in={input.r2}"
        " -S {output.r1}"
        " --read2-out {output.r2}"
        " > {log}"


rule remove_contamination:
    threads:
        8
    output:
        r1=temp("tmp/2-noadapters/{name}.1.fastq.gz"),
        r2=temp("tmp/2-noadapters/{name}.2.fastq.gz"),
    input:
        r1="tmp/1-noumi/{name}.1.fastq.gz",
        r2="tmp/1-noumi/{name}.2.fastq.gz",
    log:
        "log/2-noadapters/{name}.trimmed.log"
    group: "clean_reads"
    shell:
        "cutadapt"
        " -j {threads}"
        " -Z"
        " -e 0.15"
        " -A TTTTTCTTTTCTTTTTTCTTTTCCTTCCTTCTAA"
        " --discard-trimmed"
        " -o {output.r1}"
        " -p {output.r2}"
        " {input.r1}"
        " {input.r2}"
        " > {log}"


rule barcodes:
    """File with list of barcodes needed for demultiplexing"""
    output:
        barcodes_fasta=temp("tmp/3-barcodes/{fastqbase}.fasta")
    run:
        with open(output.barcodes_fasta, "w") as f:
            for library in libraries:
                if library.fastqbase != wildcards.fastqbase:
                    continue
                f.write(f">{library.name}\n^{library.barcode}\n")


for fastq_base, libs in fastq_map.items():

    rule:
        output:
            expand("final/fastq/{library.name}_R{read}.fastq.gz", library=libs, read=(1, 2)),
        input:
            r1="tmp/2-noadapters/{fastqbase}.1.fastq.gz".format(fastqbase=fastq_base),
            r2="tmp/2-noadapters/{fastqbase}.2.fastq.gz".format(fastqbase=fastq_base),
            barcodes_fasta="tmp/3-barcodes/{fastqbase}.fasta".format(fastqbase=fastq_base),
        params:
            r1=lambda wildcards: "final/fastq/{name}_R1.fastq.gz",
            r2=lambda wildcards: "final/fastq/{name}_R2.fastq.gz",
            fastqbase=fastq_base,
        threads: 8
        log:
            "log/3-demultiplexed/{fastqbase}.log".format(fastqbase=fastq_base)
        shell:
            "cutadapt"
            " -j {threads}"
            " -e 1"
            " --compression-level=4"
            " -g file:{input.barcodes_fasta}"
            " -o {params.r1}"
            " -p {params.r2}"
            " --discard-untrimmed"
            " {input.r1}"
            " {input.r2}"
            " > {log}"


def set_demultiplex_rule_names():
    """
    This sets the names of the demultiplexing rules, which need to be
    defined anonymously because they are defined (above) in a loop.
    """
    prefix = "tmp/2-noadapters/"
    for rul in workflow.rules:
        if not "barcodes_fasta" in rul.input.keys():
            # Ensure we get the demultiplexing rules only
            continue
        input = rul.input["r1"]
        assert input.startswith(prefix)
        # Remove the prefix and the ".1.fastq.gz" suffix
        rul.name = "demultiplex_" + rul.input["r1"][len(prefix):-11]


set_demultiplex_rule_names()


rule bowtie2:
    threads:
        19  # One fewer than available to allow other jobs to run in parallel
    output:
        bam=temp("tmp/4-mapped/{sample}_rep{replicate}.bam")
    input:
        r1="final/fastq/{sample}_rep{replicate}_R1.fastq.gz",
        r2="final/fastq/{sample}_rep{replicate}_R2.fastq.gz",
    log:
        "log/4-mapped/{sample}_rep{replicate}.log"
    # TODO
    # - --sensitive (instead of --fast) would be default
    # - write uncompressed BAM?
    # - filter unmapped reads directly? (samtools view -F 4 or bowtie2 --no-unal)
    # - add RG header
    shell:
        "bowtie2"
        " -p {threads}"
        " -x {config[bowtie_index_name]}"
        " -1 {input.r1}"
        " -2 {input.r2}"
        " --fast"
        " 2> {log}"
        " "
        "| samtools sort -@ {threads} -o {output.bam} -"


rule pool_replicates:
    output:
        bam=temp("tmp/4-mapped/{sample}_pooled.bam")
    input:
        bam_replicates=lambda wildcards: expand(
            "tmp/4-mapped/{{sample}}_rep{replicates}.bam",
            replicates=get_replicates(libraries, wildcards.sample))
    run:
        if len(input.bam_replicates) == 1:
            os.link(input.bam_replicates[0], output.bam)
        else:
            # samtools merge output is already sorted
            shell("samtools merge {output.bam} {input.bam_replicates}")


rule convert_to_single_end:
    """Convert SAM files to single-end for marking duplicates"""
    output:
        bam=temp("tmp/5-mapped_se/{library}.bam")
    input:
        bam="tmp/4-mapped/{library}.bam"
    group: "duplicate_marking"
    run:
        se_bam.convert_paired_end_to_single_end_bam(
            input.bam,
            output.bam,
            keep_unmapped=False)


# TODO have a look at UMI-tools also
rule mark_duplicates:
    """UMI-aware duplicate marking with je suite"""
    output:
        bam=temp("tmp/6-dupmarked/{library}.bam"),
        metrics="stats/6-dupmarked/{library}.metrics"
    input:
        bam="tmp/5-mapped_se/{library}.bam"
    group: "duplicate_marking"
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


rule mark_pe_duplicates:
    """Select duplicate-flagged alignments and mark them in the PE file"""
    output:
        bam=temp("tmp/7-dedup/{library}.bam")
    input:
        target_bam="tmp/4-mapped/{library}.bam",
        proxy_bam="tmp/6-dupmarked/{library}.bam"
    group: "duplicate_marking"
    run:
        se_bam.mark_duplicates_by_proxy_bam(
            input.target_bam,
            input.proxy_bam,
            output.bam)


rule remove_exclude_regions:
    output:
        bam="final/bam/{library}.bam"
    input:
        bam="tmp/7-dedup/{library}.bam",
        bed=config["exclude_regions"]
    shell:
        "bedtools"
        " intersect"
        " -v"
        " -abam {input.bam}"
        " -b {input.bed}"
        " > {output.bam}"


rule insert_size_metrics:
    output:
        txt="stats/final/{name}.insertsizes.txt",
        pdf="stats/final/{name}.insertsizes.pdf",
    input:
        bam="final/bam/{name}.bam"
    shell:
        "picard"
        " CollectInsertSizeMetrics"
        " I={input.bam}"
        " O={output.txt}"
        " HISTOGRAM_FILE={output.pdf}"
        " MINIMUM_PCT=0.5"
        " STOP_AFTER=10000000"


rule unscaled_bigwig:
    output:
        bw="final/bigwig/{library}.unscaled.bw"
    input:
        bam="final/bam/{library}.bam",
        bai="final/bam/{library}.bai",
        genome_size="stats/genome_size.txt",
    log:
        "log/final/{library}.unscaled.bw.log"
    threads: 19
    shell:
        "bamCoverage"
        " -p {threads}"
        " --normalizeUsing RPGC"
        " --effectiveGenomeSize $(< {input.genome_size})"
        " -b {input.bam}"
        " -o {output.bw}"
        " --binSize 1"
        " 2> {log}"


for group in scaling_groups:
    rule:
        output:
            factors=temp(["stats/8-factors/{library.name}.factor.txt".format(library=np.treatment) for np in group.normalization_pairs]),
            info="stats/8-scalinginfo/{scaling_group_name}.txt".format(scaling_group_name=group.name)
        input:
            treatments=["stats/final/{library.name}.flagstat.txt".format(library=np.treatment) for np in group.normalization_pairs],
            controls=["stats/final/{library.name}.flagstat.txt".format(library=np.control) for np in group.normalization_pairs],
            genome_size="stats/genome_size.txt",
        params:
            scaling_group=group,
        run:
            with open(output.info, "w") as outf:
                factors = compute_scaling(
                    params.scaling_group,
                    input.treatments,
                    input.controls,
                    outf,
                    genome_size=read_int_from_file(input.genome_size),
                    fragment_size=config["fragment_size"],
                )
                for factor, factor_path in zip(factors, output.factors):
                    with open(factor_path, "w") as f:
                        print(factor, file=f)


def set_compute_scaling_rule_names():
    """
    This sets the names of the compute scaling rules, which need to be
    defined anonymously because they are defined (above) in a loop.
    """
    prefix = "stats/8-scalinginfo/"
    for rul in workflow.rules:
        if "controls" in rul.input.keys():
            rul.name = "compute_scaling_factors_group_" + rul.output["info"][len(prefix):-4]


set_compute_scaling_rule_names()


rule summarize_scaling_factors:
    output:
        info="reports/scalinginfo.txt"
    input:
        factors=["stats/8-scalinginfo/{scaling_group_name}.txt".format(scaling_group_name=group.name) for group in scaling_groups]
    run:
        with open(output.info, "w") as f:
            print("sample_name", "#reads", "n_scaled_reads", "input_name", "n_input_reads", "factor", "scaling_group", sep="\t", file=f)
            for factor in input.factors:
                with open(factor) as factor_file:
                    for line in factor_file:
                        print(line.rstrip(), file=f)


rule extract_fragment_size:
    input:
        insertsizes="{base}.insertsizes.txt"
    output:
        fragsize=temp("{base}.fragsize.txt")
    run:
        with open(output.fragsize, "w") as f:
            print(int(parse_insert_size_metrics(input.insertsizes)["median_insert_size"]),
                  file=f)


rule scaled_bigwig:
    output:
        bw="final/bigwig/{library}.scaled.bw"
    input:
        factor="stats/8-factors/{library}.factor.txt",
        fragsize="stats/final/{library}.fragsize.txt",
        bam="final/bam/{library}.bam",
        bai="final/bam/{library}.bai",
    threads: 19
    log:
        "log/final/{library}.scaled.bw.log"
    shell:
        # TODO also run this
        # - with "--binSize 50 --smoothLength 150"
        # - with "--binSize 500 --smoothLength 5000"
        "bamCoverage"
        " -p {threads}"
        " --binSize 1"
        " --extendReads $(< {input.fragsize})"
        " --scaleFactor $(< {input.factor})"
        " --bam {input.bam}"
        " -o {output.bw}"
        " 2> {log}"


rule stats:
    output:
        txt="stats/9-stats/{library}.txt"
    input:
        mapped_flagstat="stats/4-mapped/{library}.flagstat.txt",
        metrics="stats/6-dupmarked/{library}.metrics",
        dedup_flagstat="stats/7-dedup/{library}.flagstat.txt",
        final_flagstat="stats/final/{library}.flagstat.txt",
        insertsizes="stats/final/{library}.insertsizes.txt",
    run:
        row = []
        for flagstat, name in [
            (input.mapped_flagstat, "mapped"),
            (input.dedup_flagstat, "dedup"),
            (input.final_flagstat, "final"),
        ]:
            mapped_reads = flagstat_mapped_reads(flagstat)
            row.append(mapped_reads)

        row.append(parse_duplication_metrics(input.metrics)["estimated_library_size"])
        row.append(parse_duplication_metrics(input.metrics)["percent_duplication"])
        row.append(parse_insert_size_metrics(input.insertsizes)["median_insert_size"])
        with open(output.txt, "w") as f:
            print("mapped", "dedup_mapped", "final_mapped", "library_size", "percent_duplication", "insert_size", sep="\t", file=f)
            print(*row, sep="\t", file=f)


rule stats_summary:
    output:
        txt="reports/stats_summary.txt"
    input:
        expand("stats/9-stats/{library.name}.txt", library=libraries) + expand("stats/9-stats/{pool.name}.txt", pool=pools)
    run:
        stats_summaries = [parse_stats_fields(st_file) for st_file in input]

        # I am considering we want the keys to be in a specific order
        header = [
            "library",
            "mapped",
            "dedup_mapped",
            "final_mapped",
            "library_size",
            "percent_duplication",
            "insert_size",
        ]

        with open(output.txt, "w") as f:
            print(*header, sep="\t", file=f)
            for stats_file in input:
                summary = parse_stats_fields(stats_file)
                row = [summary[k] for k in header]
                print(*row, sep="\t", file=f)


rule compute_effective_genome_size:
    output:
        txt="stats/genome_size.txt"
    input:
        fasta=config["reference_fasta"]
    run:
        with open(output.txt, "w") as f:
            print(compute_genome_size(input.fasta), file=f)


rule samtools_index:
    output:
        "{name}.bai"
    input:
        "{name}.bam"
    shell:
        "samtools index {input} {output}"


rule samtools_idxstats:
    output:
        txt="stats/{name}.idxstats.txt"
    input:
        bam="tmp/{name}.bam",
        bai="tmp/{name}.bai",
    shell:
        "samtools idxstats {input.bam} > {output.txt}"


rule samtools_flagstat:
    output:
        txt="stats/{name}.flagstat.txt"
    input:
        bam="tmp/{name}.bam"
    shell:
        "samtools flagstat {input.bam} > {output.txt}"


rule samtools_flagstat_final:
    output:
        txt="stats/final/{name}.flagstat.txt"
    input:
        bam="final/bam/{name}.bam"
    shell:
        "samtools flagstat {input.bam} > {output.txt}"
