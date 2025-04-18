#!/usr/bin/env python3

import argparse

"""
Filter unfixed mutations

Take all entries that meet the following criteria:

  1. frequency is below MAX_FREQUENCY
  2. real_read ratio is above REAL_READS_RATIO_THRESH
  3. wildtype and mutation allele min reads both above MIN_READS_THRESH
"""
DESCRIPTION = """mix_extract.py - extract unfixed mutations from forup file"""

MAX_FREQUENCY = 95
REAL_READS_RATIO_THRESH = 0.8
MIN_READS_THRESH = 4

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description=DESCRIPTION)
        parser.add_argument("forup", help="forup file")
        parser.add_argument("--max_frequency", type=float, default=MAX_FREQUENCY, help="filter frequency")
        parser.add_argument("--real_reads_ratio_thresh", type=float, default=REAL_READS_RATIO_THRESH, help="real-reads ratio threshold")
        parser.add_argument("--min_reads_thresh", type=int, default=MIN_READS_THRESH, help="min reads threshold")

        args = parser.parse_args()

        with open(args.forup) as forup:
            for line in forup:
                line = line.strip()
                a = line.split("\t")
                mut_freq = float(a[4].replace("%", ""))
                if mut_freq <= args.max_frequency:
                    b = a[7].split("=")
                    c = a[6].split("=")
                    real = float(b[0]) + float(c[0])
                    # 1. real-read ratio >= threshold
                    # 2. min reads number >= threshold for both wildtype and mutant alleles
                    if (real / float(a[5]) > args.real_reads_ratio_thresh and
                        int(b[0]) >= args.min_reads_thresh and int(c[0]) >= args.min_reads_thresh):
                            d = b[1].split(":")
                            e = c[1].split(":")
                            if int(d[0]) >= 1 and int(d[1]) >= 1 and int(e[0]) >= 1 and int(e[1]) >= 1:
                                print(line)
