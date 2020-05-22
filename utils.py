import os
import pysam
from typing import NamedTuple
from itertools import groupby, islice


class ParseError(Exception):
    pass


class Library(NamedTuple):
    sample: str
    replicate: str
    barcode: str
    fastqbase: str

    @property
    def name(self):
        return f"{self.sample}_replicate{self.replicate}"


class TreatmentControlPair(NamedTuple):
    treatment: Library
    control: Library


def read_libraries():
    for row in read_tsv("experiment.tsv"):
        yield Library(*row)


def read_controls(libraries):
    library_map = {
        (library.sample, library.replicate): library for library in libraries}
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

    # PICARD issue: it leaves library size blank sometimes, I think for small
    # sample sizes. 
    if len(values) == (len(header)-1):
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
    with open(stats_file) as f:
        header = f.readline().strip().split("\t")
        values = f.readline().strip().split("\t")
        result = {key.lower(): value for key, value in zip(header, values)}
        result["library"] = os.path.splitext(os.path.basename(stats_file))[0]
        return result

def parse_scaling_factor(filename):
    with open(filename) as f:
        return f.readline().strip()