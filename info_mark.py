#!/usr/bin/env python3

import argparse, os, math
from collections import Counter

"""info_mark.py - transform .mixfor file to .mixmark file
"""
DESCRIPTION = """info_mark.py - transform mixfor to mixmark file
"""

def chi_squared(a, b, c, d):
    if (b + d) == 0:
        return 0
    else:
        n = a + b + c + d
        return ((n * (a * d - b * c) ** 2) /
                ((a + b) * (c + d) * (a + c) * (b + d)))


if __name__ == "__main__":
        parser = argparse.ArgumentParser(description=DESCRIPTION)
        parser.add_argument("mixfor", help="mixfor file")
        args = parser.parse_args()

        k = 0
        readlen = 0

        with open(args.mixfor) as infile:
            for line in infile:
                a = line.strip().split("\t")
                a13_len = len(a[13])
                b = a[14].split(",")
                for b_elem in b:
                    i = int(b_elem)
                    if i > readlen:
                        readlen = i

                location = []
                location_str = []
                loca = 0
                dis_sum = 0
                num = 0
                med = "%.0f" % (a13_len / 2.0)  # default value

                for j in range(a13_len):
                    nu = a[13][j]
                    if nu in {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}:
                        b_j = int(b[j])
                        location_str.append(b[j])
                        location.append(b_j)
                        loca += b_j
                        #med = "%.0f" % (a13_len / 2.0)
                        dis = math.fabs(j - float(med))
                        dis_sum += dis
                        num += 1

                if num == 0:
                    aveloca = "%.0f" % 0.0
                    avedis = "%.0f" % 0.0
                else:
                    aveloca = "%.0f" % (loca / num)
                    avedis = "%.0f" % (dis_sum / num)
                ratioloca = "%.2f" % (float(aveloca) / float(readlen))
                if float(med) == 0:
                    ratiodis = "%.2f" % 0
                else:
                    ratiodis =  "%.2f" % (float(avedis) / float(med))
                c1 = a[6].strip().split("=")
                c2 = list(map(int, c1[1].split(":")))
                d1 = a[7].strip().split("=")
                d2 = list(map(int, d1[1].split(":")))
                p = "NA"

                #print("C2: ", c2)
                if (c2[0] + c2[1] != 0 and d2[0] + d2[1] != 0 and
                    c2[0] + d2[0] != 0 and c2[1] + d2[1] != 0):
                    chi = chi_squared(c2[0], c2[1], d2[0], d2[1])
                    if chi <= 2.71:
                        p = ">=0.1"
                    elif chi > 2.71 and chi <= 3.84:
                        p = "0.05~0.1"
                    elif chi > 3.84 and chi <= 6.63:
                        p = "0.01~0.05"
                    elif chi > 6.63 and chi <= 7.88:
                        p = "0.005~0.01"
                    elif chi > 7.88:
                        p = "0.005"
                #print("p: ", p)

                loc_counter = Counter()
                for m in location:
                    loc_counter[m] += 1
                # only sum the values that are > 1
                #print(loc_counter)
                all = sum([i for i in loc_counter.values() if i > 1])
                per = "%.2f" % (all / float(d1[0]))
                #print("per = ", per)
                outrow = [
                    ratioloca, "%d:%s" % (readlen, aveloca), ratiodis,
                    "%s:%s" % (med, avedis), p, per,
                    "%s:%s" % (all, d1[0]),
                    line.strip(),
                    " ".join(location_str)
                ]
                print("\t".join(outrow))
