"""
Run the Minute pipeline

Calls Snakemake to produce all the output files.

Any arguments that this wrapper script does not recognize are forwarded to Snakemake.
This can be used to provide file(s) to create, targets to run or any other Snakemake
options. For example, this runs the "full" target (including fingerprinting) in dry-run mode:

    minute run --dryrun full

Run 'snakemake --help' or see the Snakemake documentation to see valid snakemake arguments.
"""
import importlib
import logging
import re
import subprocess
import sys
from pathlib import Path
from ruamel.yaml import YAML
from ruamel.yaml import scanner

from .. import libraries_unused_in_groups, read_libraries, read_scaling_groups, make_references

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "--cores",
        "-c",
        metavar="N",
        type=int,
        help="Run on at most N CPU cores in parallel. Default: Use as many cores as available)",
    )

    parser.add_argument(
        "--dryrun",
        "-n",
        default=False,
        action="store_true",
        help="Do not execute anything",
    )


def main(args, arguments):
    run_snakemake(**vars(args), arguments=arguments)


def run_snakemake(
    dryrun=False,
    cores=None,
    arguments=None,
):
    required = [
        "references",
        "umi_length",
        "fragment_size",
        "max_barcode_errors",
    ]

    validate_config_file(Path("minute.yaml"), required)

    try:
        libraries = list(read_libraries("libraries.tsv"))
        scaling_groups = list(read_scaling_groups("groups.tsv", libraries))
    except FileNotFoundError as e:
        sys.exit(
            f"Samples configuration file '{e.filename}' not found. "
            f"Please see the documentation for how to create it."
        )

    check_fastq_basenames_exist(libraries)
    with importlib.resources.path("minute", "Snakefile") as snakefile:
        command = [
            "snakemake", f"--cores={'all' if cores is None else cores}", "-p", "-s", snakefile
        ]
        if dryrun:
            command += ["--dryrun"]
        if arguments:
            command += arguments
        logger.debug("Running: %s", " ".join(str(c) for c in command))
        exit_code = subprocess.call(command)

    warn_about_unused_libraries(libraries, scaling_groups)
    sys.exit(exit_code)


def validate_config_file(yaml, required):
    """
    Checks that the configuration minute.yaml file exists and has valid
    structure and field values.

    Arguments:
        yaml: Path to config file
        required: List of required fields
    """
    try:
        config = YAML(typ="safe").load(yaml)
    except FileNotFoundError as e:
        sys.exit(
            f"Pipeline configuration file '{e.filename}' not found. "
            f"Please see the documentation for how to create it."
        )
    except scanner.ScannerError as e:
        sys.exit(
            f"Pipeline configuration file '{yaml}' not correctly formed. "
            f"See error below.\n\n"
            f"{e}"
            f"\n\nCheck example minute.yaml file on the "
            f"documentation for reference. "
        )
    check_required_fields_exist(yaml, required)

    for name, ref in config["references"].items():
        fasta = Path(ref["fasta"])
        if not fasta.exists():
            sys.exit(
                f"Reference file {fasta} for genome {name} not found."
            )

    make_references(config["references"], config.get("aligner", "bowtie2"))

def check_required_fields_exist(yaml, required):
    """
    Checks that required fields for minute execution that have no defaults 
    are present in the YAML file. If any required field is missing it prints
    an error and exits.

    Arguments:
        yaml: Path to YAML config file
        required: List of required fields
    """
    config = YAML(typ="safe").load(yaml)
    yaml_fields = set(config.keys())
    required = set(required)
    if not required.issubset(yaml_fields):
        sys.exit(
            f"Missing required fields in '{yaml}': "
            f"{', '.join([f for f in required.difference(yaml_fields)])}.\n"
            f"Present: {', '.join([f for f in yaml_fields])}."
        )

def check_fastq_basenames_exist(libraries, fastq_dir = "fastq"):
    """
    Check that every fastq base specified in libraries.tsv is present as
    a pair of FASTQ files _R1.fastq.gz, _R2.fastq.gz

    Arguments:
        libraries: List of Library objects
        fastq_dir: Directory where FASTQ files are
    """
    basenames = set([lib.fastqbase for lib in libraries])
    missing = []
    for b in basenames:
        r1 = Path(Path(fastq_dir) / (b + "_R1.fastq.gz"))
        if not r1.exists():
            missing.append(str(r1.name))

        r2 = Path(Path(fastq_dir) / (b + "_R2.fastq.gz"))
        if not r2.exists():
            missing.append(str(r2.name))

    if len(missing) > 0:
        sys.exit(
            f"Missing files in {fastq_dir} directory with fastq base in 'libraries.tsv': "
            f"{', '.join([m for m in missing])}"
        )

def warn_about_unused_libraries(libraries, scaling_groups, limit=4):
    """
    Print warnings for libraries unused in groups.tsv.

    Arguments:
        limit: How many unused libraries to show at most
    """
    unused_libs = libraries_unused_in_groups(libraries, scaling_groups)
    if unused_libs:
        n = len(unused_libs)
        logger.warning(
            "%s",
            f"{n} librar{'y' if n == 1 else 'ies'} present in libraries.tsv "
            f"are not used in groups.tsv:"
        )
        if n == limit + 1:
            # avoid printing "... and 1 more" although
            # we could have just printed the omitted lib
            limit = n
        for unused_lib in unused_libs[:limit]:
            logger.warning("%s", f"- Library {unused_lib.sample}, replicate {unused_lib.replicate}")
        if n > limit:
            logger.warning("... and %d more", n - limit)
