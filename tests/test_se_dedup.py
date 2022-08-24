from pathlib import Path

from minute.se_dedup import AlignmentIndex


def test_index_alignment():
    alignment_file = Path("../testdata") / "mapped.bam"
    test_index = AlignmentIndex(alignment_file, 6, 5, 20)

    assert(len(test_index.r1_reads_by_position) == 1)
    assert(len(test_index.r1_reads_by_position["chr12"]) == 7)
    assert(len(test_index.r2_only_reads_by_position) == 2)
    assert(len(test_index.r2_only_reads_by_position["chr13"]) == 1)
