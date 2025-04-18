#!/usr/bin/env python3

"""
Python port of RESR_PIPELINE mix_pileup_merge.pl
"""
import argparse

DESCRIPTION = """mix_pileup_merge.py - combine unfixed file with SAM pileup"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('forfile', help=".for file")
    parser.add_argument('pileup', help="pileup file from samtools mpileup")
    args = parser.parse_args()

    pos_exists = set()
    pos2info = {}

    with open(args.forfile) as forfile:
        for line in forfile:
            line = line.strip()
            comps = line.split('\t')
            pos = int(comps[1])
            pos_exists.add(pos)
            pos2info[pos] = line

    with open(args.pileup) as pileup:
        for line in pileup:
            line = line.strip()
            b = line.split('\t')
            pos = int(b[1])
            if pos in pos_exists:
                try:
                    outrow = [pos2info[pos], b[4], b[6], b[5]]
                    print('\t'.join(outrow))
                except IndexError:
                    print("INDEX ERROR AT POS: ", pos, " b: ", b)
                    raise
