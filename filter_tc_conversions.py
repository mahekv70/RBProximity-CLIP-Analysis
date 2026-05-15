#!/usr/bin/env python3
"""
Filter fPAR-CLIP BAM alignments for diagnostic T-to-C conversions.

This script retains only mapped reads containing:
    - T>C mismatches on the forward strand
    - A>G mismatches on the reverse strand

These strand-specific substitutions correspond to crosslink-induced
T-to-C conversions after accounting for reverse-complement alignment.

Usage:
    python filter_tc_conversions.py input.bam output.bam

Requirements:
    pysam
"""

import sys
import pysam


def parse_md_tag(md_tag):
    """
    Parse the MD tag from a BAM alignment.

    The MD tag encodes matches, mismatches and deletions relative to the
    reference genome without requiring access to the reference FASTA.

    Returns:
        A list of operations:
            ("M", length)   matched reference bases
            ("X", base)     mismatched reference base
            ("D", sequence) deleted reference sequence
    """
    operations = []
    number_buffer = ""
    i = 0

    while i < len(md_tag):
        char = md_tag[i]

        if char.isdigit():
            number_buffer += char
            i += 1
            continue

        if number_buffer:
            operations.append(("M", int(number_buffer)))
            number_buffer = ""

        if char == "^":
            i += 1
            deleted_sequence = ""

            while i < len(md_tag) and md_tag[i].isalpha():
                deleted_sequence += md_tag[i]
                i += 1

            operations.append(("D", deleted_sequence))
            continue

        if char.isalpha():
            operations.append(("X", char))
            i += 1
            continue

        i += 1

    if number_buffer:
        operations.append(("M", int(number_buffer)))

    return operations


def read_has_tc_conversion(read):
    """
    Determine whether a read contains a strand-aware T>C conversion.

    For reads aligned to the forward strand, a crosslink-induced conversion
    appears as T in the reference and C in the read.

    For reads aligned to the reverse strand, the equivalent event appears
    as A in the reference and G in the read.
    """
    try:
        md_tag = read.get_tag("MD")
    except KeyError:
        return False

    sequence = read.query_sequence

    if sequence is None:
        return False

    query_position = 0

    for operation, value in parse_md_tag(md_tag):

        if operation == "M":
            query_position += value

        elif operation == "D":
            # Deletions consume reference bases but not query/read bases.
            continue

        elif operation == "X":
            if query_position >= len(sequence):
                return False

            reference_base = value.upper()
            read_base = sequence[query_position].upper()

            if not read.is_reverse:
                if reference_base == "T" and read_base == "C":
                    return True
            else:
                if reference_base == "A" and read_base == "G":
                    return True

            query_position += 1

    return False


def filter_bam_for_tc_conversions(input_bam, output_bam):
    """
    Write reads containing diagnostic T>C conversions to a new BAM file.
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out:

            for read in bam_in:
                if read.is_unmapped:
                    continue

                if read_has_tc_conversion(read):
                    bam_out.write(read)


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python filter_tc_conversions.py input.bam output.bam"
        )

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]

    filter_bam_for_tc_conversions(input_bam, output_bam)


if __name__ == "__main__":
    main()
