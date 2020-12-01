import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from io import StringIO
from itertools import groupby
from typing import List, Iterable, Dict, Optional, Tuple

from xopen import xopen


class ParseError(Exception):
    pass


@dataclass
class Reference:
    name: str
    fasta: Path
    bowtie_index: Path
    exclude_bed: Path


@dataclass
class Library:
    sample: str


@dataclass
class Replicate(Library):
    replicate: str
    barcode: str
    fastqbase: str

    @property
    def name(self):
        return f"{self.sample}_rep{self.replicate}"


@dataclass
class Pool(Library):
    replicates: List[Replicate]

    @property
    def name(self):
        return f"{self.sample}_pooled"


@dataclass
class LibraryWithReference:
    library: Library
    reference: str

    @property
    def name(self):
        return f"{self.library.name}.{self.reference}"


@dataclass
class TreatmentControlPair:
    treatment: LibraryWithReference
    control: LibraryWithReference

    @property
    def reference(self) -> str:
        assert self.treatment.reference == self.control.reference
        return self.treatment.reference


@dataclass
class ScalingGroup:
    normalization_pairs: List[TreatmentControlPair]
    name: str


def read_libraries() -> Iterable[Replicate]:
    for row in read_tsv("libraries.tsv", columns=4):
        yield Replicate(*row)


def make_pools(libraries) -> Iterable[Pool]:
    samples = defaultdict(list)
    for library in libraries:
        samples[library.sample].append(library)
    for sample, replicates in samples.items():
        yield Pool(sample=sample, replicates=replicates)


def read_scaling_groups(replicates: List[Replicate]) -> Iterable[ScalingGroup]:
    library_map: Dict[Tuple[str, str], Library] = {
        (rep.sample, rep.replicate): rep for rep in replicates
    }
    for pool in make_pools(replicates):
        library_map[(pool.sample, "pooled")] = pool

    scaling_map = defaultdict(list)
    for row in read_tsv("groups.tsv", columns=5):
        treatment_name, replicate_id, control_name, scaling_group, reference = row
        treatment_lib = library_map[(treatment_name, replicate_id)]
        control_lib = library_map[(control_name, replicate_id)]
        treatment = LibraryWithReference(treatment_lib, reference)
        control = LibraryWithReference(control_lib, reference)
        scaling_map[scaling_group].append(TreatmentControlPair(treatment, control))

    for name, normalization_pairs in scaling_map.items():
        yield ScalingGroup(normalization_pairs, name)


def flatten_scaling_groups(groups: Iterable[ScalingGroup], controls: bool = True) -> Iterable[LibraryWithReference]:
    for group in groups:
        for pair in group.normalization_pairs:
            yield pair.treatment
            if controls:
                yield pair.control


def make_references(config) -> Dict[str, Reference]:
    references = dict()
    for name, ref in config.items():
        fasta = Path(ref["fasta"])
        exclude_bed = Path(ref["exclude"])
        try:
            bowtie_index = detect_bowtie_index_name(ref["fasta"])
        except FileNotFoundError as e:
            sys.exit(str(e))
        references[name] = Reference(
            name=name,
            fasta=fasta,
            bowtie_index=bowtie_index,
            exclude_bed=exclude_bed,
        )
    return references


def read_tsv(path, columns: int) -> Iterable[List[str]]:
    """
    Read a tab-separated value file from path, allowing "#"-prefixed comments

    Yield a list of fields for every row (ignoring comments and empty lines)

    If the number of fields in a row does not match *columns*, a ParseError
    is raised.
    """
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            fields = line.strip().split("\t")
            if len(fields) != columns:
                raise ParseError(
                    f"Expected {columns} tab-separated fields in {path}, but found {len(fields)}")
            yield fields


@dataclass
class Flagstat:
    mapped_reads: int


def parse_flagstat(path) -> Flagstat:
    """Read "samtools flagstat" output and return the number of mapped reads"""
    count = None
    with open(path) as f:
        for line in f:
            if " mapped (" in line:
                count = int(line.split(maxsplit=1)[0])
                break
    return Flagstat(mapped_reads=count)


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


def compute_scaling(scaling_group, treatments, controls, infofile, genome_sizes, fragment_size):
    first = True
    scaling_factor = -1
    for pair, treatment_path, control_path, genome_size in zip(scaling_group.normalization_pairs, treatments, controls, genome_sizes):
        treatment_reads = parse_flagstat(treatment_path).mapped_reads
        control_reads = parse_flagstat(control_path).mapped_reads
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
        print(pair.treatment.name, treatment_reads, scaled_treatment_reads, pair.control.name, control_reads, sample_scaling_factor, scaling_group.name, sep="\t", file=infofile)

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
    return result


def read_int_from_file(path) -> int:
    with open(path) as f:
        data = f.read()
    return int(data.strip())


def compute_genome_size(fasta: str) -> int:
    n = 0
    with xopen(fasta) as f:
        for line in f:
            if line[:1] != ">":
                line = line.rstrip().upper()
                n += len(line) - line.count("N")
    return n


def detect_bowtie_index_name(fasta_path: str) -> Path:
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
        f"No Bowtie2 index found for '{fasta_path}', expected one of\n- "
        + "\n- ".join(str(b) + bowtie_index_extension for b in bases)
    )


def get_replicates(libraries, sample):
    replicates = [lib.replicate for lib in libraries if lib.sample == sample]
    return replicates


def get_normalization_pairs(scaling_groups) -> List[TreatmentControlPair]:
    return [pair for group in scaling_groups for pair in group.normalization_pairs]


def format_metadata_overview(references, libraries, maplibs, scaling_groups) -> str:
    f = StringIO()

    print("# References", file=f)
    for name, reference in references.items():
        print(" -", reference, file=f)
    print(file=f)

    print("# Libraries", file=f)
    for library in libraries:
        print(" -", library, file=f)

    print(file=f)
    print("# Pools", file=f)
    for maplib in maplibs:
        if isinstance(maplib.library, Pool):
            pool = maplib.library
            print(" -", pool.name, "(replicates:", ", ".join(r.replicate for r in pool.replicates) + ")", file=f)

    print(file=f)
    print("# Scaling groups", file=f)

    for group in scaling_groups:
        print("# Group", group.name, "- Normalization Pairs (treatment -- control)", file=f)
        for pair in group.normalization_pairs:
            print(" -", pair.treatment.library.name, "--", pair.control.library.name,
                "(reference: {})".format(pair.reference), file=f)
    return f.getvalue()


def is_snakemake_calling_itself() -> bool:
    return "snakemake/__main__.py" in sys.argv[0]


def map_fastq_prefix_to_list_of_libraries(replicates: List[Replicate]) -> Dict[str, List[Replicate]]:
    return {
        fastq_base: list(reps)
        for fastq_base, reps in
        groupby(sorted(replicates, key=lambda rep: rep.fastqbase), key=lambda rep: rep.fastqbase)
    }
