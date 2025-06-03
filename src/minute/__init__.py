import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from io import StringIO
from itertools import groupby
from typing import List, Iterable, Dict, Tuple, Optional, Union, Set

from xopen import xopen

try:
    from importlib.metadata import version as _version
except ImportError:
    from importlib_metadata import version as _version

__version__ = _version("minute")


class ParseError(Exception):
    pass


@dataclass
class Reference:
    name: str
    fasta: Path
    bowtie_index: Path
    exclude_bed: Optional[Path]


@dataclass(eq=True, frozen=True)
class Library:
    sample: str


@dataclass(eq=True, frozen=True)
class Replicate(Library):
    replicate: str
    fastqbase: str

    @property
    def name(self):
        return f"{self.sample}_rep{self.replicate}"

    def has_umi(self) -> bool:
        with xopen(f"fastq/{self.fastqbase}_R1.fastq.gz") as f:
            line = f.readline()
        return bool(re.search("_[ACGTNacgtn]+", line))


@dataclass(eq=True, frozen=True)
class MultiplexedReplicate(Replicate):
    barcode: str

    def has_umi(self) -> bool:
        # Demultiplexing and UMI removal are done at the same time in the pipeline, so
        # we know that a replicate that needed demultiplexing has UMIs
        return True


@dataclass(eq=True, frozen=True)
class Pool(Library):
    replicates: List[Replicate]

    @property
    def name(self):
        return f"{self.sample}_pooled"

    def has_umi(self) -> bool:
        return all(lib.has_umi() for lib in self.replicates)


@dataclass(eq=True, frozen=True)
class LibraryWithReference:
    library: Library
    reference: str

    @property
    def name(self):
        return f"{self.library.name}.{self.reference}"


@dataclass(eq=True, frozen=True)
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


def read_libraries(path: Union[os.PathLike, str]) -> Iterable[Replicate]:
    for row in read_tsv(path, columns=4):
        if row[2] == ".":
            yield Replicate(sample=row[0], replicate=row[1], fastqbase=row[3])
        else:
            yield MultiplexedReplicate(
                sample=row[0], replicate=row[1], barcode=row[2], fastqbase=row[3]
            )


def make_pools(libraries) -> Iterable[Pool]:
    samples = defaultdict(list)
    for library in libraries:
        samples[library.sample].append(library)
    for sample, replicates in samples.items():
        yield Pool(sample=sample, replicates=replicates)


def read_scaling_groups(
        path: Union[os.PathLike, str], replicates: List[Replicate]
) -> Iterable[ScalingGroup]:

    library_map: Dict[Tuple[str, str], Library] = {
        (rep.sample, rep.replicate): rep for rep in replicates
    }
    for pool in make_pools(replicates):
        library_map[(pool.sample, "pooled")] = pool

    scaling_map = defaultdict(list)
    for row in read_tsv(path, columns=5):
        pairs = []
        treatment_name, replicate_id, control_name, scaling_group, reference = row
        treatment_lib = library_map[(treatment_name, replicate_id)]
        control_lib = library_map[(control_name, replicate_id)]
        pairs.append((treatment_lib, control_lib))

        if isinstance(treatment_lib, Pool):
            assert isinstance(control_lib, Pool)
            assert len(treatment_lib.replicates) == len(control_lib.replicates)
            for tlib, clib in zip(treatment_lib.replicates, control_lib.replicates):
                pairs.append((tlib, clib))

        for treatment_lib, control_lib in pairs:
            treatment = LibraryWithReference(treatment_lib, reference)
            control = LibraryWithReference(control_lib, reference)
            pair = TreatmentControlPair(treatment, control)
            if pair not in scaling_map.get(scaling_group, []):
                scaling_map[scaling_group].append(pair)

    for name, normalization_pairs in scaling_map.items():
        yield ScalingGroup(normalization_pairs, name)


def flatten_scaling_groups(groups: Iterable[ScalingGroup], controls: bool = True) -> Iterable[LibraryWithReference]:
    seen = set()

    for group in groups:
        for pair in group.normalization_pairs:
            maplibs = [pair.treatment]
            if controls:
                maplibs.append(pair.control)
            for maplib in maplibs:
                if maplib.name not in seen:
                    seen.add(maplib.name)
                    yield maplib


def get_all_controls(groups: Iterable[ScalingGroup]) -> Iterable[LibraryWithReference]:
    seen = set()
    maplibs = list()
    for group in groups:
        for pair in group.normalization_pairs:
            if pair.control.name not in seen:
                seen.add(pair.control.name)
                yield pair.control


def get_all_pools(maplibs: Iterable[LibraryWithReference]) -> List[LibraryWithReference]:
    return [m for m in maplibs if isinstance(m.library, Pool)]


def get_all_replicates(maplibs: Iterable[LibraryWithReference]) -> List[LibraryWithReference]:
    return [m for m in maplibs if not isinstance(m.library, Pool)]

def get_maplib_by_name(maplibs: Iterable[LibraryWithReference], name: str) -> LibraryWithReference:
    for m in maplibs:
        if m.library.name == name:
            return m

def make_references(config, aligner) -> Dict[str, Reference]:
    references = dict()
    for name, ref in config.items():
        fasta = Path(ref["fasta"])
        if not fasta.exists():
            sys.exit(f"Reference file {fasta} not found.")
        exclude_bed = Path(ref["exclude"]) if ref["exclude"] else None
        bowtie_index = None
        if aligner == "bowtie2":
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
    Read a whitespace-separated value file from path, allowing "#"-prefixed comments

    Yield a list of fields for every row (ignoring comments and empty lines)

    If the number of fields in a row does not match *columns*, a ParseError
    is raised.
    """
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            fields = line.strip().split()
            if len(fields) != columns:
                raise ParseError(
                    f"Expected {columns} whitespace-separated fields in {path}, "
                    f"but found {len(fields)}. Including a space character in a value "
                    f"can cause this error."
                )
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

    # Je outputs missing value as ?
    if result.get("percent_duplication") == "?":
        result["percent_duplication"] = "NA"

    return result


def compute_scaling(scaling_group, treatments, controls, infofile, genome_sizes, fragment_size):
    first = True
    scaling_factor = -1
    for pair, treatment_path, control_path, genome_size in zip(scaling_group.normalization_pairs, treatments, controls, genome_sizes):
        treatment_reads = parse_flagstat(treatment_path).mapped_reads
        control_reads = parse_flagstat(control_path).mapped_reads

        try:
            if first:
                scaling_factor = (
                    genome_size / fragment_size / treatment_reads * control_reads
                )
                treatment_reads_ref = treatment_reads
                control_reads_ref = control_reads
                first = False

            sample_scaling_factor = scaling_factor / control_reads
            scaled_treatment_reads = sample_scaling_factor * treatment_reads
        except ZeroDivisionError:
            sample_scaling_factor = "NA"
            print(pair.treatment.name, treatment_reads, scaled_treatment_reads, pair.control.name, control_reads, sample_scaling_factor, scaling_group.name, sep="\t", file=infofile)
            yield sample_scaling_factor
        else:
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


def get_sample_replicates(libraries, sample):
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


def map_fastq_prefix_to_list_of_libraries(
        replicates: List[MultiplexedReplicate]
) -> Dict[str, List[MultiplexedReplicate]]:
    return {
        fastq_base: list(reps)
        for fastq_base, reps in
        groupby(sorted(replicates, key=lambda rep: rep.fastqbase), key=lambda rep: rep.fastqbase)
    }


def libraries_unused_in_groups(libraries: List[Replicate], groups: List[ScalingGroup]) -> List[Library]:
    unused: Set[Replicate] = set(libraries)
    for group in groups:
        for pair in group.normalization_pairs:
            for lib in pair.treatment.library, pair.control.library:
                if isinstance(lib, Replicate):
                    unused.discard(lib)
                elif isinstance(lib, Pool):
                    unused = {a for a in unused if a.sample != lib.sample}
                else:
                    assert False, "Expected only Replicate or Pool"
    return list(unused)


def lander_waterman(x, c, n):
    """
    Lander-Waterman equation states:
    C/X = 1 - exp( -N/X )
    where
     X = number of distinct molecules in library
     N = number of read pairs
     C = number of distinct fragments observed in read pairs

    Returns y value of the function: y = c/x -1 + exp(-n/x)
    """
    return c/x - 1 + math.exp(-n / x)


def estimate_library_size(total_reads, duplicate_reads):
    """
    Picard procedure for estimating library size translated directly to Python.
    See Picard repository for the original: 
    https://github.com/broadinstitute/picard/blob/5295289523f8526b42a08b6e0f0111c8ed5e4399/src/main/java/picard/sam/DuplicationMetrics.java#L143
    """
    unique_reads = total_reads - duplicate_reads
    if total_reads == 0 or duplicate_reads == 0:
        return "NA"

    m = 1.0
    M = 100.0

    if (unique_reads >= total_reads or lander_waterman(m * unique_reads, unique_reads, total_reads) < 0):
        raise ValueError(f"Invalid values for pairs and unique pairs: {total_reads} {unique_reads}")

    # find value of M, large enough to act as other side for bisection method
    while lander_waterman(M * unique_reads, unique_reads, total_reads) > 0:
        M *= 10.0

    for i in range(40):
        r = (m + M) / 2.0
        u = lander_waterman(r*unique_reads, unique_reads, total_reads)

        if u == 0:
            break
        elif u > 0:
            m = r
        elif u < 0:
            M = r

    return unique_reads * (m + M) / 2.0
