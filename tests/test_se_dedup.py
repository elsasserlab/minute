import os
from pathlib import Path

from minute.se_dedup import ContigIndex, FirstMateDeduplicator, PysamAlignmentFile


def test_index_alignment():
    alignment_file = Path("testdata") / "mapped.bam"
    alignment_index = Path("testdata") / "mapped.bai"
    alignment = PysamAlignmentFile(alignment_file, alignment_index)
    test_index = ContigIndex(alignment, "chr12", 6, 5)

    assert(len(test_index.r1_reads_by_position) == 7)
    assert(len(test_index.r2_only_reads_by_position) == 3)


def test_se_dedup(tmp_path):
    alignment_file = Path("testdata") / "mapped.bam"
    alignment_index = Path("testdata") / "mapped.bai"
    dedup_file = tmp_path / "dedup.bam"
    deduplicator = FirstMateDeduplicator(umi_length=6,
                                         multimap_cutoff=5,
                                         stub_length=20)

    dups_summary = deduplicator.deduplicate(alignment_file, alignment_index, dedup_file)

    assert(dups_summary.total == 12)
    assert(dups_summary.total_dups == 1)
    assert(dups_summary.r1_dups == 1)
    assert(dups_summary.r2_only_dups == 0)
    assert(os.path.exists(dedup_file))

    os.unlink(dedup_file)
