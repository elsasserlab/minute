#!/usr/bin/env python3
"""
Run the Minute pipeline

Calls Snakemake to produce all the output files.

Any arguments that this wrapper script does not recognize are forwarded to Snakemake.
This can be used to provide file(s) to create, targets to run or any other Snakemake
options. For example, this runs the "quick" target (without fingerprinting) in dry-run mode:

    minute run --dryrun quick

Run 'snakemake --help' or see the Snakemake documentation to see valid snakemake arguments.
"""
import logging
import subprocess
import sys
from argparse import ArgumentParser
from pathlib import Path

from ruamel.yaml import YAML
import importlib.resources

logger = logging.getLogger(__name__)


def main(arguments=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = ArgumentParser(description=__doc__, prog="minute")
    parser.add_argument("--debug", action="store_true", help="Print some debugging information")
    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("run", help="Run the pipeline")

    arg = subparser.add_argument
    subparser.set_defaults(func=run_snakemake)
    arg(
        "--cores",
        "-c",
        metavar="N",
        type=int,
        help="Run on at most N CPU cores in parallel. Default: Use as many cores as available)",
    )

    arg(
        "--dryrun",
        "-n",
        default=False,
        action="store_true",
        help="Do not execute anything",
    )
    args, remainder = parser.parse_known_args(arguments)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    if hasattr(args, "func"):
        subcommand = args.func
        del args.func
    else:
        parser.error("Please provide the name of a subcommand to run")
    subcommand(**vars(args), arguments=remainder)


def run_snakemake(
    dryrun=False,
    cores=None,
    arguments=None,
):
    try:
        _ = YAML(typ="safe").load(Path("config.yaml"))
    except FileNotFoundError as e:
        sys.exit(
            f"Pipeline configuration file '{e.filename}' not found. "
            f"Please see the documentation for how to create it."
        )
    with importlib.resources.path("minute", "Snakefile") as snakefile:
        command = [
            "snakemake", f"--cores={'all' if cores is None else cores}", "-p", "-s", snakefile
        ]
        if dryrun:
            command += ["--dryrun"]
        if arguments:
            command += arguments
        logger.debug("Running: %s", " ".join(str(c) for c in command))
        sys.exit(subprocess.call(command))


if __name__ == "__main__":
    main()
