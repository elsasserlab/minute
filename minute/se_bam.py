import argparse
import os
import pysam
import tempfile


def _is_read_valid(alignment, keep_unmapped):
    return alignment.is_read1 and (keep_unmapped or not alignment.is_unmapped)


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

    # Silence spurious "Could not retrieve index file" warnings
    verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(source_bam) as source_file,\
         pysam.AlignmentFile(out, 'wb', header=source_file.header) as out_file:
        pysam.set_verbosity(verbosity)
        for alignment in source_file:
            if alignment.query_name in dup_ids:
                alignment.is_duplicate = True

            if not filter_dups or not alignment.is_duplicate:
                out_file.write(alignment)


def _get_duplicate_ids(bam):
    # Silence spurious "Could not retrieve index file" warnings
    verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'rb') as bam_file:
        pysam.set_verbosity(verbosity)
        idlist = {alignment.query_name
            for alignment in bam_file if alignment.is_duplicate}
        return idlist
