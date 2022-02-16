#!/usr/bin/env python3
"""
Run the Minute pipeline

Calls Snakemake to produce all the output files.
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

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("run", help="Run the pipeline")

    arg = subparser.add_argument
    subparser.set_defaults(func=run_snakemake)
    arg(
        "--dryrun",
        "-n",
        default=False,
        action="store_true",
        help="Do not execute anything",
    )
    arg(
        "--cores",
        "--jobs",
        "-j",
        metavar="N",
        type=int,
        default=0,
        help="Run on at most N CPU cores in parallel. Default: Use as many cores as available)",
    )
    arg(
        "targets",
        nargs="*",
        default=[],
        help="File(s) to create or targets to run. If omitted, the full pipeline is run.",
    )
    args = parser.parse_args(arguments)
    if hasattr(args, "func"):
        subcommand = args.func
        del args.func
    else:
        parser.error("Please provide the name of a subcommand to run")
    subcommand(**vars(args))


def run_snakemake(
    dryrun=False,
    cores=0,
    targets=None,
):
    try:
        config = YAML(typ="safe").load(Path("config.yaml"))
        print("config", config)
    except FileNotFoundError as e:
        sys.exit(
            f"Pipeline configuration file '{e.filename}' not found. "
            f"Please see the documentation for how to create it.")
    with importlib.resources.path("minute", "Snakefile") as snakefile:
        command = ["snakemake", f"--cores={'all' if cores == 0 else cores}", "-p", "-s", snakefile]
        if dryrun:
            command += ["--dryrun"]
        if targets:
            command += targets
        sys.exit(subprocess.call(command))


if __name__ == "__main__":
    main()
