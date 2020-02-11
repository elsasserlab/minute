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
    with pysam.AlignmentFile(bam, 'rb') as infile,\
         pysam.AlignmentFile(out, 'wb', header=infile.header) as outfile:
        for alignment in infile:
            if is_read_valid(alignment, keep_unmapped):
                new_alignment = convert_to_single_end(alignment)
                outfile.write(new_alignment)

def get_duplicate_ids(bam):
    idlist = []
    with pysam.AlignmentFile(bam, 'rb') as bamfile:
        for alignment in bamfile:
            if alignment.is_duplicate:
                idlist.append(alignment.id)
    return idlist

def mark_duplicates_by_proxy_bam(target_bam, proxy_bam, out, filter_dups=True):
    """
    Iterates through a pair of bam files, assumed to have the same
    number of reads, sorted the same way. target_bam reads are marked as dups
    if they are marked as such in the proxy_bam file.

    Keyword argument:
    filter_dups -- Remove duplicates in the output file.
    """
    with pysam.AlignmentFile(target_bam) as target_file,\
         pysam.AlignmentFile(proxy_bam) as proxy_file,\
         pysam.AlignmentFile(out, 'rb', header=target_bam.header) as outfile:

        for proxy_alignment in proxy_file:
            target_alignment = target_file.next()
            target_alignment.is_duplicate = proxy_alignment.is_duplicate

            if not target_alignment.is_duplicate or not filter_dups:
                outfile.write(target_alignment)

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
