import os
import pysam
from typing import NamedTuple
from itertools import groupby, islice
from pathlib import Path


class ParseError(Exception):
    pass


class Library(NamedTuple):
    sample: str
    replicate: str
    barcode: str
    fastqbase: str

    @property
    def name(self):
        rep_str = f"replicate{self.replicate}"
        if self.replicate == "pooled":
            rep_str = f"pooled"

        return f"{self.sample}_{rep_str}"


class TreatmentControlPair(NamedTuple):
    treatment: Library
    control: Library


def read_libraries():
    for row in read_tsv("experiment.tsv"):
        yield Library(*row)


def read_controls(libraries):
    library_map = {
        (library.sample, library.replicate): library for library in libraries}

    for library in libraries:
        library_map[(library.sample, 'pooled')] = Library(
            sample=library.sample,
            replicate='pooled',
            barcode='-',
            fastqbase='-')

    for row in read_tsv("controls.tsv"):
        treatment = library_map[(row[0], row[1])]
        control = library_map[(row[2], row[1])]
        yield TreatmentControlPair(treatment, control)


def read_tsv(path):
    """
    Read a tab-separated value file from path, allowing "#"-prefixed comments

    Yield a list of fields for every row (ignoring comments and empty lines)
    """
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            fields = line.strip().split("\t")
            yield fields


def flagstat_mapped_reads(path):
    """Read "samtools flagstat" output and return the number of mapped reads"""
    with open(path) as f:
        for line in f:
            if " mapped (" in line:
                return int(line.split(maxsplit=1)[0])


def parse_insert_size_metrics(path):
    return parse_picard_metrics(path, metrics_class="picard.analysis.InsertSizeMetrics")


def parse_duplication_metrics(path):
    return parse_picard_metrics(path, metrics_class="picard.sam.DuplicationMetrics")


def parse_picard_metrics(path, metrics_class: str):
    """
    Parse the METRICS section in a metrics file created by Picard. The string
    given after 'METRICS CLASS' in the file must match the metrics_class
    parameter. Anything else in the file is not parsed.

    Return a dictionary that has keys such as "median_insert_size",
    "estimated_library_size" etc. (this depends on the headers in the file).
    """
    def float_or_int(s):
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return s

    with open(path) as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                fields = line.strip().split(sep="\t")
                if fields[1] != metrics_class:
                    raise ParseError(
                        "While parsing Picard metrics file '{}':"
                        "Expected metrics class {}, but found {}".format(
                            path, metrics_class, fields[1]))
                break
        header = next(f).strip().split("\t")
        values = next(f).strip().split("\t")

    # Picard issue: it leaves library size blank sometimes, I think for small
    # sample sizes. 
    if len(values) == len(header) - 1:
        values.append('NA')

    result = {key.lower(): float_or_int(value) for key, value in zip(header, values)}
    return result


def compute_scaling(normalization_pairs, treatments, controls, infofile, genome_size, fragment_size):
    print("sample_name", "#reads", "n_scaled_reads", "input_name", "n_input_reads", "factor", sep="\t", file=infofile)
    first = True
    scaling_factor = -1
    for pair, treatment_path, control_path in zip(normalization_pairs, treatments, controls):
        treatment_reads = flagstat_mapped_reads(treatment_path)
        control_reads = flagstat_mapped_reads(control_path)
        if first:
            scaling_factor = (
                genome_size / fragment_size / treatment_reads * control_reads
            )
            treatment_reads_ref = treatment_reads
            control_reads_ref = control_reads
            first = False

        sample_scaling_factor = scaling_factor / control_reads
        scaled_treatment_reads = sample_scaling_factor * treatment_reads

        # TODO factor this out
        print(pair.treatment.name, treatment_reads, scaled_treatment_reads, pair.control.name, control_reads, sample_scaling_factor, sep="\t", file=infofile)

        # TODO scaled.idxstats.txt file

        yield sample_scaling_factor


def parse_stats_fields(stats_file):
    """
    Parse contents of a stats file created in the stats rule, which consists
    of a header line and a values line.

    Return a dictionary where keys are defined by the values of the header.
    """
    with open(stats_file) as f:
        header = f.readline().strip().split("\t")
        values = f.readline().strip().split("\t")
        result = {key.lower(): value for key, value in zip(header, values)}
        result["library"] = Path(stats_file).stem
        return result


def detect_bowtie_index_name(fasta_path):
    """
    Given the path to a reference FASTA (which may optionally be compressed),
    detect the base name of the Bowtie2 index assumed to be in the same
    location.

    Given "ref.fasta.gz", this function checks for the existence of
    - "ref.fasta.gz.1.bt2"
    - "ref.fasta.1.bt2"
    - "ref.1.bt2"
    in that order and then returns the name of the first file it finds
    minus the ".1.bt2" suffix.
    """
    path = Path(fasta_path)
    bowtie_index_extension = ".1.bt2"
    bases = [path]
    if path.suffix == ".gz":
        bases.append(path.with_suffix(""))
    if bases[-1].suffix:
        bases.append(bases[-1].with_suffix(""))
    for base in bases:
        if base.with_name(base.name + bowtie_index_extension).exists():
            return base
    raise FileNotFoundError(
        "No Bowtie2 index found, expected one of\n- "
        + "\n- ".join(str(b) + bowtie_index_extension for b in bases)
    )


def get_replicates(libraries, sample):
    replicates = [lib.replicate for lib in libraries if lib.sample==sample]
    return replicates


def infer_pooled_libraries(libraries):
    samples = set([library.sample for library in libraries])
    return([
        Library(sample=sample, replicate='pooled', barcode='', fastqbase='')
        for sample in samples
    ])
