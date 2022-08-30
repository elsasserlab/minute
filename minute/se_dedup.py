import datetime
import math
import os
import pysam
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterator, List
from umi_tools import UMIClusterer


def count_sequences(seq_list):
    counts = defaultdict(int)
    for s in seq_list:
        counts[s] += 1
    return counts


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


class PysamAlignmentFile:
    def __init__(self, bam_file, index_file):
        self.bam = bam_file
        self.index = index_file

    def get_mapped_segments(self, contig: str = None, min_quality: int = 0, max_quality: int = 44) -> Iterator[pysam.AlignedSegment]:
        verbosity = pysam.set_verbosity(0)
        with pysam.AlignmentFile(self.bam, mode='rb',
                                 index_filename=self.index) as handle:
            pysam.set_verbosity(verbosity)
            for read in handle.fetch(contig=contig):
                if not read.is_unmapped and read.mapping_quality >= min_quality and read.mapping_quality < max_quality:
                    yield read

    def get_mapped_count(self) -> int:
        """Returns total number of fragments in the whole alignment"""
        # Discordant alignments can be counted twice if done per contig.
        seen = []
        with pysam.AlignmentFile(self.bam, mode='rb',
                                 index_filename=self.index) as handle:
            for read in handle:
                if not read.is_unmapped:
                    seen.append(read.query_name)

        return len(set(seen))

    def get_alignment_statistics(self) -> List:
        with pysam.AlignmentFile(self.bam, mode='rb', index_filename=self.index) as handle:
            return handle.get_index_statistics()

    def get_contig_names(self) -> Iterator[str]:
        with pysam.AlignmentFile(self.bam, mode='rb', index_filename=self.index) as handle:
            for name in handle.references:
                yield name

    def write(self, dst_file, exclude_ids, umi_length=6):
        verbosity = pysam.set_verbosity(0)
        with pysam.AlignmentFile(self.bam, mode='rb') as src, \
                pysam.AlignmentFile(dst_file, 'wb', header=src.header) as dst:
            pysam.set_verbosity(verbosity)
            for r in src:
                if not r.query_name[:-umi_length - 1] in exclude_ids:
                    dst.write(r)


class SegmentSummary:
    name: str
    umi: bytes
    stub: bytes

    def __init__(self, read_info: pysam.AlignedSegment, umi_length: int, stub_length: int):
        self.name = self.parse_name(read_info, umi_length)
        self.umi = self.parse_umi(read_info, umi_length)
        self.stub = str.encode(read_info.query_sequence[0:stub_length])
        self.mapping_quality = read_info.mapping_quality
        self.score = self.get_score(read_info)

    @staticmethod
    def parse_name(read, umi_length):
        return read.query_name[:-umi_length - 1]

    @staticmethod
    def parse_umi(read, umi_length):
        return str.encode(read.query_name[-umi_length:])

    @staticmethod
    def get_score(read):
        return sum([int(q) for q in read.query_alignment_qualities])

    def is_multi(self, multimap_cutoff):
        if self.mapping_quality < multimap_cutoff:
            return True


class AlignedSegmentSummary(SegmentSummary):
    def __init__(self, read_info: pysam.AlignedSegment, umi_length: int, stub_length: int):
        super().__init__(read_info, umi_length, stub_length)
        self.mate = self.get_mate(read_info)
        self.start = read_info.reference_start
        self.mate_start = read_info.next_reference_start
        # We don't need the stub for aligned reads
        self.stub = b""

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


class ContigIndex:
    """Indexes a contig BAM alignment in two categories:
    - R1 reads: Mapped reads that have mapping quality >= multimap_cutoff
    and are flagged as R1. These can be part of a proper pair or not, since
    the deduplication is done based solely on R1 sequence. They are indexed
    by genomic position.
    - R2-only reads: Mapped reads that have mapping quality >= multimap_cutoff,
    are flagged as R2 and whose matching R1 is not found in the alignment, or
    is unmapped. They are indexed by genomic position.
    """
    def __init__(self, alignment: PysamAlignmentFile, contig: str, umi_length: int, multimap_cutoff: int):
        self.alignment = alignment
        self.r1_reads_by_position = defaultdict(list)
        self.r2_only_reads_by_position = defaultdict(list)

        self.umi_length = umi_length
        self.multimap_cutoff = multimap_cutoff

        seen_1st = self._index_first_mates(contig)
        seen_2nd = self._index_second_mates(contig, seen_1st)

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

    def _index_first_mates(self, contig):
        seen = list()
        for read in self.alignment.get_mapped_segments(contig, self.multimap_cutoff):
            aln = AlignedSegmentSummary(read,
                                        umi_length=self.umi_length,
                                        stub_length=0)
            if aln.mate == 1:
                if not aln.is_multi(multimap_cutoff=self.multimap_cutoff):
                    pos = read.reference_start
                    self.r1_reads_by_position[pos].append(aln)
                    seen.append(aln.name)
        return set(seen)

    def _index_second_mates(self, contig, seen):
        second_ids = []
        for read in self.alignment.get_mapped_segments(contig, self.multimap_cutoff):
            aln = AlignedSegmentSummary(read, umi_length=self.umi_length,
                                        stub_length=0)

            if aln.mate == 2 and aln.name not in seen:
                if not aln.is_multi(multimap_cutoff=self.multimap_cutoff):
                    pos = read.reference_start
                    self.r2_only_reads_by_position[pos].append(aln)
                    second_ids.append(aln.name)
        return set(second_ids)


@dataclass
class DedupSummary:
    library: str
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
        header = "## METRICS CLASS\tminute.DuplicationMetrics"
        fields = ["LIBRARY",
                  "READ_PAIRS_EXAMINED",
                  "READ_PAIR_DUPLICATES",
                  "R1_DUPLICATES",
                  "R2_ONLY_DUPLICATES",
                  "MULTIMAPPING_DUPLICATES",
                  "PERCENT_DUPLICATION",
                  "ESTIMATED_LIBRARY_SIZE"]
        values = [self.library,
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


def print_timed(s):
    time_str = str(datetime.datetime.now())
    print(f"{time_str}: {s}")


class FirstMateDeduplicator:
    """Deduplicates a BAM file using R1 position only when possible.

    Read pairs where only R2 is mapped are deduplicated using their position.
    Multimap pers are all pooled together and clustered hierarchically,
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

    def deduplicate_contig(self, alignment, contig):
        index = ContigIndex(alignment, contig, self.umi_length, self.multimap_cutoff)
        r1_dups = self._find_duplicate_ids(index.r1_reads_by_position)
        r2_dups = self._find_duplicate_ids(index.r2_only_reads_by_position)
        return [r1_dups, r2_dups]

    def deduplicate(self, src_bam, scr_index, dst_bam):
        alignment = PysamAlignmentFile(src_bam, scr_index)
        all_mapped_dups = dict()
        for contig in alignment.get_contig_names():
            print_timed(f"Dedup contig {contig}")
            all_mapped_dups[contig] = self.deduplicate_contig(alignment, contig)

        # Deduplicate multimappers
        multi = [SegmentSummary(r, umi_length=self.umi_length, stub_length=self.stub_length)
                 for r in alignment.get_mapped_segments(max_quality = self.multimap_cutoff)]
        print_timed(f"Deduplicating {len(multi)} multimapping reads")
        multi_dups = self.mark_duplicates_by_umi_and_sequence(multi)

        print_timed("Writing without the duplicates")
        alignment.write(dst_bam, all_mapped_dups, self.umi_length)
        print_timed("Done")
        contig_dups = self.summarise_contig_dups(alignment, all_mapped_dups)
        total = alignment.get_mapped_count()

        return DedupSummary(
            library=os.path.basename(alignment.bam),
            total=total,
            total_dups=sum(contig_dups)+len(multi_dups),
            r1_dups=contig_dups[0],
            r2_only_dups=contig_dups[1],
            multi_dups=len(multi_dups)
        )

    @staticmethod
    def summarise_contig_dups(alignment, dups_dict):
        r1_dups = 0
        r2_dups = 0

        for contig in dups_dict:
            r1_dups += len(dups_dict[contig][0])
            r2_dups += len(dups_dict[contig][1])

        return [r1_dups, r2_dups]

    def _find_duplicate_ids(self, aln_index):
        duplicates = list()
        for i, dup_candidates in aln_index.items():
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
