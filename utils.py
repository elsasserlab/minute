from typing import NamedTuple
from itertools import groupby


class Library(NamedTuple):
    sample: str
    replicate: str
    barcode: str
    fastqbase: str

    @property
    def name(self):
        return f"{self.sample}_replicate{self.replicate}"


def read_experiment_description():
    with open("experiment.tsv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            fields[1] = fields[1].lstrip("R")
            yield Library(*fields)

