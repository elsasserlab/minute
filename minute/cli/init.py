"""
Create and initialize a new pipeline directory
"""
import logging
import os
from pathlib import Path
import importlib.resources
from typing import Optional

from . import CommandLineError


logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "--reads",
        type=Path,
        help="Raw reads directory with paired-end FASTQ files."
    )
    parser.add_argument("directory", type=Path, help="New pipeline directory to create")


def main(args, arguments):
    if arguments:
        raise CommandLineError("These arguments are unknown: %s", arguments)
    run_init(**vars(args))


def run_init(directory: Path, reads: Optional[Path]):
    if " " in str(directory):
        raise CommandLineError("The name of the pipeline directory must not contain spaces")

    if reads is not None and not reads.is_dir():
        raise CommandLineError(f"'{reads}' must be a directory")

    try:
        directory.mkdir()
    except OSError as e:
        raise CommandLineError(e)

    if reads is not None:
        relative_symlink(reads, directory / "fastq")
    else:
        logger.info(
            "Option --reads not used, please create and populate directory %s/fastq/ manually",
            directory,
        )

    configuration = importlib.resources.read_text("minute", "minute.yaml")
    with open(Path(directory) / "minute.yaml", "w") as f:
        f.write(configuration)

    logger.info("Pipeline directory %s created", directory)
    logger.info(
        'Edit %s/%s and run "cd %s && minute run" to start the analysis',
        directory,
        "minute.yaml",
        directory,
    )


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

