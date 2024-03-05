"""
Create and initialize a new pipeline directory
"""
import logging
import os
import shutil
from pathlib import Path
import importlib.resources
from typing import Optional

from . import CommandLineError

from minute import read_tsv

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "--reads",
        type=Path,
        help="Raw reads directory with paired-end FASTQ files."
    )
    parser.add_argument(
        "--barcodes",
        type=Path,
        help="Barcodes description file."
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Name of the input sample (FASTQ name without _R1/_R2)."
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional minute.yaml file to use as configuration. If not provided, a template will be created."
    )
    parser.add_argument("directory", type=Path, help="New pipeline directory to create")


def main(args, arguments):
    if arguments:
        raise CommandLineError("These arguments are unknown: %s", arguments)
    run_init(**vars(args))


def run_init(directory: Path, reads: Optional[Path], barcodes: Optional[Path], input: Optional[str], config: Optional[Path]):
    if " " in str(directory):
        raise CommandLineError("The name of the pipeline directory must not contain spaces")

    if reads is not None and not reads.is_dir():
        raise CommandLineError(f"'{reads}' must be a directory")

    if reads is None and barcodes is not None:
        raise CommandLineError("--reads parameter must be specified if --barcodes option is used")
        logger.info(
            "Option --reads not used, please create and populate directory %s/fastq/ manually",
            directory,
        )

    if barcodes is not None:
        if not os.path.isfile(barcodes):
            raise CommandLineError(f"--barcodes '{barcodes}' file not found")
        if input is None:
            raise CommandLineError("--input must be specified if --barcodes option is used")

        libraries = make_libraries_from_barcodes_and_reads(barcodes, reads)
        try:
            groups = make_groups_from_barcodes_and_reads(barcodes, reads, input)
        except ValueError as e:
            raise CommandLineError(f"Invalid --input value: {e}")

    if config is not None:
        if not os.path.isfile(config):
            raise CommandLineError(f"--config '{config}' file not found")
        

    try:
        directory.mkdir()
    except OSError as e:
        raise CommandLineError(e)

    if reads is not None:
        relative_symlink(reads, directory / "fastq")

    if barcodes is not None:
        write_tsv(libraries, directory / "libraries.tsv")
        write_tsv(groups, directory / "groups.tsv")

    if config is not None:
        shutil.copyfile(config, directory / "minute.yaml")
    else:
        configuration = importlib.resources.read_text("minute", "minute.yaml")
        with open(Path(directory) / "minute.yaml", "w") as f:
            f.write(configuration)

    logger.info("Pipeline directory %s created", directory)
    logger.info(
        'Edit %s/%s if necessary and run "cd %s && minute run" to start the analysis',
        directory,
        "minute.yaml",
        directory,
    )


def parse_library_prefixes(reads_directory):
    return [r.name[:-12] for r in reads_directory.glob("*_R1.fastq.gz")]


def make_libraries_from_barcodes_and_reads(barcodes, reads_directory):
    """
    Infers libraries info based on the fastq directory and a barcodes info file.
    Each condition in the barcodes file will be applied to each fastq.gz pair
    in reads_directory.
    """
    prefixes = parse_library_prefixes(reads_directory)
    barcodes_info = list(read_tsv(barcodes, columns=4))

    libraries_info = []
    for p in prefixes:
        for bc in barcodes_info:
            libraries_info.append(
                [f"{p}_{bc[0]}", bc[1], bc[2], p]
            )

    return libraries_info


def make_groups_from_barcodes_and_reads(barcodes, reads_directory, input):
    """
    Infers groups info based on the fastq directory, a barcodes info file and
    the name of the input sample making the following assumptions: every 
    FASTQ R1/R2 pair contains a scaling group; the reference condition used
    is the corresponding to the first line in the barcodes file, pooling the
    samples; all barcodes in the barcodes info file are applied to each
    FASTQ R1/R2 pair in the reads_directory; input must correspond to one file
    pair in the reads_directory.
    """
    prefixes = parse_library_prefixes(reads_directory)
    if input not in prefixes:
        msg = f"input ({input}) must match a FASTQ file prefix, found: {prefixes}"
        raise ValueError(msg)

    prefixes = [p for p in prefixes if p != input]
    barcodes_info = list(read_tsv(barcodes, columns=4))
    barcodes_refmap = {bc[0]: bc[3] for bc in barcodes_info}

    # Reference condition is first
    ref_condition = barcodes_info[0][0]
    conditions = list(dict.fromkeys([bc[0] for bc in barcodes_info]))
    conditions = [c for c in conditions if c != ref_condition]

    groups_info = []

    for p in prefixes:
        groups_info.append([
            f"{p}_{ref_condition}",
            "pooled",
            f"{input}_{ref_condition}",
            p,
            barcodes_refmap[ref_condition]
        ])
        for c in conditions:
            groups_info.append(
                [f"{p}_{c}", "pooled", f"{input}_{c}", p, barcodes_refmap[c]]
            )

    return groups_info


def write_tsv(values, path):
    """
    Write a tab-separated file.
    """
    with open(path, "w") as f:
        for row in values:
            f.write("\t".join([str(v) for v in row]) + "\n")


def relative_symlink(src, dst, force: bool = False):
    """
    Create a symbolic link in any directory.

    force -- if True, then overwrite an existing file/symlink
    """
    if force:
        try:
            os.remove(dst)
        except FileNotFoundError:
            pass
    target = os.path.relpath(os.path.abspath(src), start=os.path.dirname(dst))
    os.symlink(target, dst)
