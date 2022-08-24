from pathlib import Path

from minute.se_dedup import AlignmentIndex, FirstMateDeduplicator


def test_index_alignment():
    alignment_file = Path("testdata") / "mapped.bam"
    test_index = AlignmentIndex(alignment_file, 6, 5, 20)

    assert(len(test_index.r1_reads_by_position) == 1)
    assert(len(test_index.r1_reads_by_position["chr12"]) == 7)
    assert(len(test_index.r2_only_reads_by_position) == 2)
    assert(len(test_index.r2_only_reads_by_position["chr13"]) == 1)


def test_se_dedup(tmp_path):
    alignment_file = Path("testdata") / "mapped.bam"
    dedup_file = tmp_path / "dedup.bam"
    deduplicator = FirstMateDeduplicator(umi_length=6,
                                         multimap_cutoff=5,
                                         stub_length=20)

    dups_summary = deduplicator.deduplicate(alignment_file)

    assert(dups_summary.total == 12)
    assert(dups_summary.total_dups == 1)
    assert(dups_summary.r1_dups == 1)
    assert(dups_summary.r2_only_dups == 0)
