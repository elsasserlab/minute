import math
import os
import pysam
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Iterator, List
from umi_tools import UMIClusterer


def count_sequences(seq_list):
    counts = defaultdict(int)
    for s in seq_list:
        counts[s] += 1
    return counts


def get_mapped_segments(bam: os.PathLike) -> Iterator[pysam.AlignedSegment]:
    """Iterates through a BAM file returning only mapped AlignedSegment
    objects """
    verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'rb') as handle:
        pysam.set_verbosity(verbosity)
        for read in handle:
            if not read.is_unmapped:
                yield read


def write_bam_excluding_reads(src_file, dst_file, exclude_ids, umi_length=6):
    verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(src_file, 'rb') as src, \
            pysam.AlignmentFile(dst_file, 'wb', header=src.header) as dst:
        pysam.set_verbosity(verbosity)
        for r in src:
            if not r.query_name[:-umi_length - 1] in exclude_ids:
                dst.write(r)


def lander_waterman(x: int, c: int, n: int) -> float:
    """Lander-Waterman equation C/X = 1 - exp( -N/X )
            y = C/X -1 + exp(-N/X)
    """
    return c / x - 1 + math.exp(-n / x)


def estimate_library_size(total, unique, iter=40) -> float|None:
    """Estimates the size of a library based on Lander Waterman equation,
    same way as PICARD does (MIT-licensed), based on the total number of
    reads and the number of unique reads observed.

    NOTE:
    For PICARD, total reads are: PAIRS EXAMINED - OPTICAL DUPS
                unique reads are: PAIRS EXAMINED - PAIR DUPS.
    We do not estimate optical duplicates, so these are included in the
    library size estimate.
    """
    if total > 0 and unique > 0:
        lower_bound = 1.0
        upper_bound = 100.0

        wt_lower = lander_waterman(lower_bound * unique, unique, total)
        if unique >= total or wt_lower < 0:
            print(f"Non valid values for pairs and unique pairs: {total} {unique}")
            return None

        # find value of upper bound
        while lander_waterman(upper_bound * unique, unique, total) > 0:
            upper_bound *= 10.0

        for i in range(iter):
            pivot = (lower_bound + upper_bound) / 2.0
            u = lander_waterman(pivot * unique, unique, total)

            if u == 0:
                break
            elif u > 0:
                lower_bound = pivot
            elif u < 0:
                upper_bound = pivot

        return unique * (lower_bound + upper_bound) / 2.0
    else:
        return None


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


@dataclass
class DedupSummary:
    total: int
    total_dups: int
    r1_dups: int
    r2_only_dups: int
    multi_dups: int

    def fraction_duplication(self) -> float:
        return self.total_dups / self.total

    def libsize(self) -> float:
        return estimate_library_size(self.total, self.total - self.total_dups)

    def format_stats(self) -> str:
        # Each element in a pair is only counted once, so theoretically
        # we are still counting "pairs"
        header = "## METRICS CLASS\tpicard.sam.DuplicationMetrics"
        fields = ["LIBRARY",
                  "READ_PAIRS_EXAMINED",
                  "READ_PAIR_DUPLICATES",
                  "R1_DUPLICATES",
                  "R2_ONLY_DUPLICATES",
                  "MULTIMAPPING_DUPLICATES",
                  "PERCENT_DUPLICATION",
                  "ESTIMATED_LIBRARY_SIZE"]
        values = ["Unknown library",
                  self.total,
                  self.total_dups,
                  self.r1_dups,
                  self.r2_only_dups,
                  self.multi_dups,
                  self.fraction_duplication(),
                  self.libsize()]
        values_str = "\t".join(str(v) for v in values)
        fields_str = "\t".join(fields)
        return f"{header}\n{fields_str}\n{values_str}"


class FirstMateDeduplicator:
    """Deduplicates a BAM file using R1 position only when possible.

    Read pairs where only R2 is mapped are deduplicated using their position.
    Multimappers are all pooled together and clustered hierarchically,
    first by UMI, then by a sequence stub of max length stub_length, after
    which reads are deduplicated.
    """
    def __init__(self, umi_length: int = 6,
                 multimap_cutoff: int = 5,
                 stub_length: int = 20,
                 umi_mismatches: int = 1,
                 seq_mismatches: int = 2):

        self.umi_length = umi_length
        self.multimap_cutoff = multimap_cutoff
        self.stub_length = stub_length
        self.umi_mismatches = umi_mismatches
        self.seq_mismatches = seq_mismatches

    def deduplicate(self, src_bam, dst_bam):
        index = AlignmentIndex(src_bam, self.umi_length, self.multimap_cutoff,
                               self.stub_length)
        r1_dups = self._find_duplicate_ids(index.r1_reads_by_position)
        r2_dups = self._find_duplicate_ids(index.r2_only_reads_by_position)
        multi_duplicates = self.mark_duplicates_by_umi_and_sequence(
            index.multimappers)
        all_dups = set(r1_dups + r2_dups + multi_duplicates)
        write_bam_excluding_reads(src_bam, dst_bam, all_dups, self.umi_length)

        return DedupSummary(
            total=index.total,
            total_dups=len(all_dups),
            r1_dups=len(r1_dups),
            r2_only_dups=len(r2_dups),
            multi_dups=len(multi_duplicates))

    def _find_duplicate_ids(self, aln_index):
        duplicates = list()
        for ref in aln_index.keys():
            for i, dup_candidates in aln_index[ref].items():
                if len(dup_candidates) > 1:
                    dup_ids = self.mark_duplicates_by_umi(
                        dup_candidates, self.umi_mismatches)
                    duplicates.extend(dup_ids)
        return duplicates

    def mark_duplicates_by_umi(self, read_list, umi_mismatches=1):
        seqs = [r.umi for r in read_list]
        clusters = self.cluster_sequences(seqs, mismatches=umi_mismatches)
        unique_reads = []
        for c in clusters:
            # First element in the list is considered the original according
            # to UMI-tools documentation
            unique_candidates = self.filter_by_attribute(read_list, c[0], "umi")
            if len(unique_candidates) > 1:
                unique_reads.append(self.pick_best(unique_candidates))
            else:
                unique_reads.append(unique_candidates[0])

        unique_reads_set = set(r.name for r in unique_reads)
        return [r.name for r in read_list if r.name not in unique_reads_set]

    def mark_duplicates_by_umi_and_sequence(self, read_list, umi_mismatches=1, seq_mismatches=2):
        # First group by similar UMIs
        # Then for each group cluster by sequences that match all those UMIs
        # (we are not picking yet)
        seq_list = [r.umi for r in read_list]
        umi_clusters = self.cluster_sequences(seq_list, mismatches=umi_mismatches)

        unique_reads = []
        for c in umi_clusters:
            cluster_reads = self.filter_by_attribute(read_list, c, "umi")
            # Require stub to reach length - UMIClusterer requires all counts
            # to have the same length
            cluster_seqs = [r.stub for r in cluster_reads if len(r.stub) == self.stub_length]
            sub_clusters = self.cluster_sequences(cluster_seqs, mismatches=seq_mismatches)
            for s in sub_clusters:
                # First stub in the cluster is considered the "original"
                # according to UMItools documentation
                sub_cluster_reads = self.filter_by_attribute(cluster_reads, s[0], "stub")
                best_read = self.pick_best(sub_cluster_reads)
                unique_reads.append(best_read)

        unique_reads_set = set(r.name for r in unique_reads)
        return [r.name for r in read_list if r.name not in unique_reads_set]

    @staticmethod
    def cluster_sequences(seq_list, mismatches=1, cluster_method="directional"):
        clusters = []
        if len(seq_list) > 0:
            clusterer = UMIClusterer(cluster_method=cluster_method)
            counts = count_sequences(seq_list)
            clusters = clusterer(counts, threshold=mismatches)
        return clusters

    @staticmethod
    def filter_by_attribute(reads: List[AlignedSegmentSummary], allowed: List[str], attribute: str) -> List[AlignedSegmentSummary]:
        return [r for r in reads if getattr(r, attribute) in allowed]

    @staticmethod
    def pick_best(readlist):
        best = readlist[0]
        for r in readlist[1:]:
            if r.score > best.score:
                best = r
        return best
