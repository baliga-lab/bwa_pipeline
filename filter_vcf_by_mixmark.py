#!/usr/bin/env python3

"""
Retain the lines from the VCF file that are in the mixmark
file.
"""


import vcfpy
import argparse


DESCRIPTION = """filter_vcf_by_mixmark.py - filter VCF by mixmark
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('vcf', help="input VCF file")
    parser.add_argument('mixmark', help="input mixmark file")
    parser.add_argument('outfile', help="output VCF file")
    args = parser.parse_args()

    # TODO: we need to annotate the reads of the ALT allele vs the number of
    # reads on the REF allele
    included = set()
    with open(args.mixmark) as infile:
        for line in infile:
            included.add(int(line.strip().split()[8]))

    reader = vcfpy.Reader.from_path(args.vcf)
    writer = vcfpy.Writer.from_path(args.outfile, reader.header)
    for record in reader:
        if record.POS in included:
            writer.write_record(record)
