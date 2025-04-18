#!/usr/bin/env python3

import argparse
import re

"""
Resr pipeline
.mix -> .mixfor file conversion
"""

DESCRIPTION = """forup_format.py - transform forup file
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("forup", help="forup file")
    args = parser.parse_args()

    k = 0
    with open(args.forup) as infile:

        for line in infile:
            line = line.strip()
            a = line.split("\t")
            # filter variants with varscan P value failed
            # note that variants from VCF files that were run with a
            # real pvalue in varscan will have different values
            # that's why we can't use the filter in the same way
            #if a[8].startswith("Pass=1E0"):
            if a[8].startswith("PASS"):
                seq = ""
                # remove '$' and '^' characters
                a[13] = a[13].replace("$", "")
                a[13] = re.sub(r"\^.", "", a[13])
                a13_len = len(a[13])
                for i in range(a13_len):
                    a13_i = a[13][i]
                    #print("Processing '%s'" % a13_i)
                    if a13_i == '+' or a13_i == '-':
                        a1 = a[13][i + 1]
                        a2 = a[13][i + 2]
                        # k = how many bases to skip
                        if a1.isdigit() and not a2.isdigit():
                            k = int(a1) + 1
                        elif a1.isdigit() and a2.isdigit():
                            k = int(a1) * 10 + int(a2) + 2
                    elif a[13][i] != '+' and a[13][i] != '-' and k > 0:
                        k -= 1
                    elif a[13][i] != '+' and a[13][i] != '-' and k == 0:
                        seq += a[13][i]

                len1 = len(seq)
                b = a[14].split(",")
                len2 = len(b)
                #print("len1: %d len2: %d" % (len1, len2))
                if len1 == len2:
                    outrow = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], seq, a[14]]
                else:
                    outrow = ["ERROR%s" % a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], seq, a[14]]
                    #print(outrow)
                print("\t".join(outrow))
