#!/usr/bin/env python3

"""
Get the SRA file URL for the sample id
"""
import argparse


DESCRIPTION = "get_sra_url.py - get SRA file URL for download"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('download_csv', help="path to download csv file")
    parser.add_argument('sample_id', help="sample id")
    args = parser.parse_args()

    download_map = {}
    with open(args.download_csv) as infile:
        for line in infile:
            sample, url = line.strip().split('\t')
            download_map[sample] = url
    try:
        print(download_map[args.sample_id])
    except:
        raise
        exit(1)
