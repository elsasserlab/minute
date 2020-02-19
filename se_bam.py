import argparse
import os
import pysam
import tempfile


def convert_paired_end_to_single_end_bam(bam, out, keep_unmapped=False):
    """
    Iterates through a paired-end BAM file and keeps only the first
    mate of each pair.

    Keyword arguments:
    keep_unmapped -- Keep unmapped reads
    """
    with pysam.AlignmentFile(bam, 'rb') as in_file,\
         pysam.AlignmentFile(out, 'wb', header=in_file.header) as out_file:
        for alignment in in_file:
            if _is_read_valid(alignment, keep_unmapped):
                new_alignment = _convert_to_single_end(alignment)
                out_file.write(new_alignment)


def _is_read_valid(alignment, keep_unmapped):
    return alignment.is_read1 and (keep_unmapped or not alignment.is_unmapped)


def _convert_to_single_end(alignment):
    alignment.is_paired = False
    alignment.is_read1 = False
    alignment.is_proper_pair = False
    alignment.mate_is_reverse = False
    alignment.mate_is_unmapped = False
    alignment.next_reference_id = 0
    alignment.next_reference_start = 0
    alignment.next_reference_name = None
    return alignment


def mark_duplicates_by_proxy_bam(source_bam, proxy_bam, out, filter_dups=True):
    """
    Iterates through a pair of BAM files. source_bam reads are marked as
    dups if they are marked as such in the proxy_bam file.

    Keyword argument:
    filter_dups -- Remove duplicates in the output file.
    """

    dup_ids = _get_duplicate_ids(proxy_bam)

    if not dup_ids:
        os.link(source_bam, out)
        return

    with pysam.AlignmentFile(source_bam) as source_file,\
         pysam.AlignmentFile(out, 'wb', header=source_file.header) as out_file:

        for alignment in source_file:
            if alignment.query_name in dup_ids:
                alignment.is_duplicate = True

            if not filter_dups or not alignment.is_duplicate:
                out_file.write(alignment)


def _get_duplicate_ids(bam):
    with pysam.AlignmentFile(bam, 'rb') as bam_file:
        idlist = {alignment.query_name
            for alignment in bam_file if alignment.is_duplicate}
        return idlist
