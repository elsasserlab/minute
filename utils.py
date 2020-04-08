import pysam
from typing import NamedTuple
from itertools import groupby, islice


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
    """
    Parse a metrics file created by Picard InsertSizeMetrics. The histogram
    in the file is not parsed.

    Return a dictionary that has keys "median_insert_size", "mode_insert_size"
    etc.
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
                break
        header = next(f).strip().split()
        values = next(f).strip().split()
    result = {key.lower(): float_or_int(value) for key, value in zip(header, values)}
    return result


def compute_scaling(treatments, controls, infofile, genome_size, fragment_size):
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

        print(pair.treatment.name, treatment_reads, scaled_treatment_reads, pair.control.name, control_reads, sample_scaling_factor, sep="\t", file=infofile)

        # TODO scaled.idxstats.txt file

        yield sample_scaling_factor
