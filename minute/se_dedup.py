import os
import pysam
from collections import defaultdict
from enum import Enum
from typing import Iterator


def get_mapped_segments(bam: os.PathLike) -> Iterator[pysam.AlignedSegment]:
    """Iterates through a BAM file returning only mapped AlignedSegment
    objects """
    with pysam.AlignmentFile(bam) as handle:
        for read in handle:
            if not read.is_unmapped:
                yield read


class AlignmentType(Enum):
    SINGLE = 0
    MULTI = 1


class AlignedSegmentSummary:
    """Information about a pysam.AlignedSegment relevant for the
    deduplication method"""
    def __init__(self,
                 read_info: pysam.AlignedSegment,
                 umi_length: int,
                 multimap_cutoff: int,
                 stub_length: int):

        self.name = self.get_name(read_info, umi_length)
        self.umi = self.get_umi(read_info, umi_length)
        self.mate = self.get_mate(read_info)
        self.reference = read_info.reference_name
        self.start = read_info.reference_start
        self.mate_start = read_info.next_reference_start
        self.mapping_quality = read_info.mapping_quality
        self.stub = self.get_stub(read_info, multimap_cutoff, stub_length)
        self.score = self.get_score(read_info)
        self.type = self.get_type(read_info, multimap_cutoff)

    @staticmethod
    def get_type(read, multimap_cutoff):
        t = AlignmentType.SINGLE
        if read.mapping_quality < multimap_cutoff:
            t = AlignmentType.MULTI
        return t

    @staticmethod
    def get_name(read, umi_length):
        return read.query_name[:-umi_length - 1]

    @staticmethod
    def get_umi(read, umi_length):
        return str.encode(read.query_name[-umi_length:])

    @staticmethod
    def get_mate(read):
        if read.is_read1:
            return 1
        elif read.is_read2:
            return 2
        else:
            # There are reads that are not assigned a mate (flags 0 and 16).
            # This might be because if bowtie2 fails to map a pair then it tries
            # to map the reads independently.
            return -1

    @staticmethod
    def get_stub(read, multimap_cutoff, length):
        if read.mapping_quality < multimap_cutoff:
            stub = str.encode(read.query_sequence[:length])
            return stub
        else:  # Save memory storing only this info if necessary (multimappers)
            return b""

    @staticmethod
    def get_score(read):
        return sum([int(q) for q in read.query_alignment_qualities])


class AlignmentIndex:
    """Indexes a BAM alignment in three categories:
    - R1 reads: Mapped reads that have mapping quality >= multimap_cutoff
    and are flagged as R1. These can be part of a proper pair or not, since
    the deduplication is done based solely on R1 sequence. They are indexed
    by genomic position.
    - R2-only reads: Mapped reads that have mapping quality >= multimap_cutoff,
    are flagged as R2 and whose matching R1 is not found in the alignment, or
    is unmapped. They are indexed by genomic position.
    - Multimappers: Mapped reads that have mapping quality < multimap_cutoff.
    R1 reads based on the same criteria as before, and R2 where R1 is not
    found are gathered together, as their deduplication is not based on
    position, only on sequence.
    """
    def __init__(self, bam, umi_length, multimap_cutoff, stub_length):
        self.r1_reads_by_position = defaultdict(lambda: defaultdict(list))
        self.r2_only_reads_by_position = defaultdict(lambda: defaultdict(list))
        self.multimappers = list()

        self.umi_length = umi_length
        self.multimap_cutoff = multimap_cutoff
        self.stub_length = stub_length

        seen_1st = self._index_first_mates(bam)
        seen_2nd = self._index_second_mates(bam, seen_1st)

        self.r1_counts = len(seen_1st)
        self.r2_counts = len(seen_2nd)
        self.total = self.r1_counts + self.r2_counts

    @staticmethod
    def count_loci(alignment_index):
        n = 0
        for contig in alignment_index.values():
            n += len(contig)
        return n

    def print_summary(self):
        r1_npos = self.count_loci(self.r1_reads_by_position)
        r2_npos = self.count_loci(self.r2_only_reads_by_position)
        r1_ref = len(self.r1_reads_by_position)
        r2_ref = len(self.r2_only_reads_by_position)
        print(f"R1: {self.r1_counts} ({r1_npos} pos, {r1_ref} ref)")
        print(f"R2: {self.r2_counts} ({r2_npos} pos, {r2_ref} ref)")
        print(f"Multi: {len(self.multimappers)}")

    def _index_first_mates(self, bam):
        seen = list()
        for read in get_mapped_segments(bam):
            aln = AlignedSegmentSummary(read,
                                        self.umi_length,
                                        self.multimap_cutoff,
                                        self.stub_length)
            if aln.mate == 1:
                if aln.type == AlignmentType.SINGLE:
                    refname = read.reference_name
                    pos = read.reference_start
                    self.r1_reads_by_position[refname][pos].append(aln)
                    seen.append(aln.name)
                elif aln.type == AlignmentType.MULTI:
                    self.multimappers.append(aln)
                    seen.append(aln.name)
        return set(seen)

    def _index_second_mates(self, bam, seen):
        second_ids = []
        for read in get_mapped_segments(bam):
            aln = AlignedSegmentSummary(read,
                                        self.umi_length,
                                        self.multimap_cutoff,
                                        self.stub_length)

            if aln.mate == 2 and aln.name not in seen:
                if aln.type == AlignmentType.SINGLE:
                    refname = read.reference_name
                    pos = read.reference_start
                    self.r2_only_reads_by_position[refname][pos].append(aln)
                    second_ids.append(aln.name)
                elif aln.type == AlignmentType.MULTI:
                    self.multimappers.append(aln)
                    second_ids.append(aln.name)
        return set(second_ids)
