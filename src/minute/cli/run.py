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
import subprocess
import sys
from pathlib import Path
from ruamel.yaml import YAML

from .. import libraries_unused_in_groups, read_libraries, read_scaling_groups

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
    try:
        _ = YAML(typ="safe").load(Path("minute.yaml"))
    except FileNotFoundError as e:
        sys.exit(
            f"Pipeline configuration file '{e.filename}' not found. "
            f"Please see the documentation for how to create it."
        )
    try:
        libraries = list(read_libraries("libraries.tsv"))
        scaling_groups = list(read_scaling_groups("groups.tsv", libraries))
    except FileNotFoundError as e:
        sys.exit(e)

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
