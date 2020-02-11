import argparse
import pysam
import sys
import tempfile

from shutil import copyfile

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
            if is_read_valid(alignment, keep_unmapped):
                new_alignment = convert_to_single_end(alignment)
                out_file.write(new_alignment)

def get_duplicate_ids(bam):
    idlist = []

    with pysam.AlignmentFile(bam, 'rb') as bam_file:
        for alignment in bam_file:
            if alignment.is_duplicate:
                idlist.append(alignment.query_name)
    return idlist

def find_next_duplicate(alignment_file, read_id, out_file):

    alignment = next(alignment_file)
    while alignment.query_name != read_id:
        # These reads are not duplicates
        out_file.write(alignment)
        alignment = next(alignment_file)
    return alignment

def mark_duplicates_by_proxy_bam(target_bam, proxy_bam, out, filter_dups=True):
    """
    Iterates through a pair of bam files, assumed to have IDs in the same order
    (they come from the same sorted bam file). target_bam reads are marked as
    dups if they are marked as such in the proxy_bam file.

    Keyword argument:
    filter_dups -- Remove duplicates in the output file.
    """

    dup_ids = get_duplicate_ids(proxy_bam)

    if not dup_ids:
        copyfile(target_bam, out)
        return

    with pysam.AlignmentFile(target_bam) as target_file,\
         pysam.AlignmentFile(out, 'wb', header=target_file.header) as out_file:

        for read_id in dup_ids:
            target_alignment = find_next_duplicate(
                target_file,
                read_id,
                out_file)

            target_alignment.is_duplicate = True
            if not filter_dups:
                out_file.write(target_alignment)

        for remaining_alignment in target_file:
            out_file.write(remaining_alignment)

def is_read_valid(alignment, keep_unmapped):
    if alignment.is_read1:
        if keep_unmapped:
            return True
        elif not alignment.is_unmapped:
            return True

    return False

def convert_to_single_end(alignment):
    alignment.is_paired = False
    alignment.is_read1 = False
    alignment.is_proper_pair = False
    alignment.mate_is_reverse = False
    alignment.mate_is_unmapped = False
    alignment.next_reference_id = 0
    alignment.next_reference_start = 0
    alignment.next_reference_name = None
    return alignment
