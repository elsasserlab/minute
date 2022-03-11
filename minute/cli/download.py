"""
Download SRA data for accessions listed in libraries.tsv

Run this in a pipeline directory containing "libraries.tsv".
The downloaded files are written to the fastq/ directory.
"""
import logging
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from .. import read_libraries, MultiplexedReplicate

from . import CommandLineError


logger = logging.getLogger(__name__)


def add_arguments(_parser):
    pass


def main(args, arguments):
    if arguments:
        raise CommandLineError("These arguments are unknown: %s", arguments)
    run_download(**vars(args))


def run_download():
    with TemporaryDirectory(dir=".") as tmpdir:
        tmp_path = Path(tmpdir)
        for library in read_libraries("libraries.tsv"):
            if isinstance(library, MultiplexedReplicate):
                continue
            accession = library.fastqbase
            fastq = Path("fastq")
            fastq.mkdir(exist_ok=True)
            r1 = fastq / f"{accession}_R1.fastq.gz"
            r2 = fastq / f"{accession}_R2.fastq.gz"
            if r1.exists() and r2.exists():
                logger.info("Files for accession %s exist, skipping", accession)
                continue
            command = [
                "fastq-dump",
                "--outdir",
                tmpdir,
                "--gzip",
                "--split-3",
                "--defline-qual",
                "+",
                "--defline-seq",
                "@$sn",
                accession,
            ]
            logger.info("Running %s", " ".join(command))
            subprocess.check_call(command)
            (tmp_path / f"{accession}_1.fastq.gz").rename(r1)
            (tmp_path / f"{accession}_2.fastq.gz").rename(r2)
