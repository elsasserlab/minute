import argparse
import pysam
import sys
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

def mark_duplicates_by_proxy_bam(target_bam, proxy_bam, out, filter_dups=True):
    """
    Iterates through a pair of bam files, assumed to have the same
    number of reads, sorted the same way. target_bam reads are marked as dups
    if they are marked as such in the proxy_bam file.

    Keyword argument:
    filter_dups -- Remove duplicates in the output file.
    """
    #TODO: Iterate in a smarter way: These bams are sorted, so it should not
    # be necessary to check the whole list for each read.
    dup_ids = get_duplicate_ids(proxy_bam)

    with pysam.AlignmentFile(target_bam) as target_file,\
         pysam.AlignmentFile(out, 'wb', header=target_file.header) as out_file:

        for target_alignment in target_file:
            if target_alignment.query_name in dup_ids:
                target_alignment.is_duplicate = True

            if not target_alignment.is_duplicate or not filter_dups:
                out_file.write(target_alignment)

def is_read_valid(alignment, keep_unmapped):
    if alignment.is_read1:
        if keep_unmapped:
            return True
        elif not alignment.is_unmapped:
            return True

    return False

def strand_flag(alignment):
    if alignment.is_reverse:
        return 16
    else:
        return 0

def convert_to_single_end(alignment):
    alignment.flag = strand_flag(alignment)
    alignment.is_paired = False
    alignment.is_read1 = False
    alignment.is_proper_pair = False
    alignment.mate_is_reverse = False
    alignment.mate_is_unmapped = False
    alignment.next_reference_id = 0
    alignment.next_reference_start = 0
    alignment.next_reference_name = None
    return alignment


if __name__ == '__main__':
    convert_paired_end_to_single_end_bam(sys.argv[1], 'test.bam')
    # mark_duplicates_by_proxy_bam(sys.argv[1], sys.argv[2], 'test.bam')
