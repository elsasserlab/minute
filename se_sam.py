import argparse
import tempfile


def convert_paired_end_sam_to_single_end(sam, out, keep_unmapped=False, sep='\t'):
    """
    Manually iterates through a paired-end SAM file and keeps only the first
    mate of each pair. Flags and mate-pair fields are changed accordingly.

    Keyword arguments:
    keep_unmapped -- Keep unmapped reads (flag 0x04)
    sep -- Separator between fields (this defaults to tab as in SAM specs)
    
    """

    fi = open(sam)
    fo = open(out, 'w')
    flags_field = 1

    for line in fi:
        if is_header(line):
            fo.write(line)

        else:
            fields = line.split(sep)
            flags = fields[flags_field]
            if valid_flag(flags, keep_unmapped):
                new_alignment_fields = convert_to_single_end(fields)
                fo.write(sep.join(new_alignment_fields))

    fo.close()
    fi.close()

def is_header(line):
    if line.startswith('@'):
        return True
    return False

def valid_flag(flags, keep_unmapped):
    # Bits for unmapped and second in pair need to be zero
    discard_mask = 0x84
    if keep_unmapped:
        discard_mask = 0x80

    try:
        flag_int = int(flags)
        if not (flag_int & discard_mask):
            return True
        else:
            return False

    except ValueError:
        msg = "Invalid flag field value found: {}".format(flags)
        raise ValueError(msg)
        # print(msg)

def mask_paired_end_flags(flag):
    if int(flag) & 0x10: # reverse
        return "16"
    else:
        return "0"

def convert_to_single_end(fields):
    flag_field = 1
    read_id = fields[0]
    new_flag = mask_paired_end_flags(fields[flag_field])

    # mate_id, mate_start, mate_end replaced with "missing" values
    mate_info = ["*", "0", "0"]

    new_sam_fields = [read_id, new_flag] + fields[2:6] + mate_info + fields[9:]
    return new_sam_fields
