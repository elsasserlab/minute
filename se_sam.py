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
    with pysam.AlignmentFile(bam, "rb") as infile:
        with pysam.AlignmentFile(out, 'wb', header=infile.header) as outfile:
            for alignment in infile:
                if is_read_valid(alignment, keep_unmapped):
                    new_alignment = convert_to_single_end(alignment)
                    outfile.write(new_alignment)

def is_header(line):
    return line.startswith('@')

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
    alignment.next_reference_id = 0
    alignment.next_reference_start = 0
    return alignment


if __name__ == '__main__':
    convert_paired_end_to_single_end_bam(sys.argv[1], 'test.bam')
