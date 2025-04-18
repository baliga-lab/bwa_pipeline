#!/usr/bin/env python3

"""
Replacement script for PE_IS_filt.py. This rewrite is to make the
intention more clear.
support VCF files, which have a different format
"""
import argparse
import pandas
import vcfpy

DESCRIPTION = """varscan_drop_excluded.py - drop excluded regions from Varscan file"""

def process_varscan_file(varscan_path, out_path, excluded):
    df = pandas.read_csv(varscan_path, sep='\t')
    result_indices = []
    for index, row in df.iterrows():
        if not row['Position'] in excluded:
            result_indices.append(index)
    result_df = df[df.index.isin(result_indices)]
    result_df.to_csv(out_path, sep='\t', index=False)


def process_vcf_file(vcf_path, out_path, excluded):
    """We won't even try to generate a varscan file here. We just return the
    filtered VCF lines"""
    reader = vcfpy.Reader.from_path(vcf_path)
    writer = vcfpy.Writer.from_path(out_path, reader.header)
    for record in reader:
        if not record.POS in excluded:
            writer.write_record(record)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('exclude_list', help="exclude list")
    parser.add_argument('varscan', help="input varscan file")
    parser.add_argument('outfile', help="output varscan file")
    args = parser.parse_args()
    excluded = set()
    with open(args.exclude_list) as infile:
        for line in infile:
            excluded.add(int(line.strip()))

    if args.varscan.endswith('.vcf'):
        process_vcf_file(args.varscan, args.outfile, excluded)
    else:
        process_varscan_file(args.varscan, args.outfile, excluded)
