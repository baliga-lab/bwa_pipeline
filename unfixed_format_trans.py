#!/usr/bin/env python3

"""
Transform Varscan format to unfixed format
"""
import argparse
import vcfpy

DESCRIPTION = """unfixed_format_trans.py - transform Varscan to unfixed format"""


def process_vcf_file(vcf_path):
    reader = vcfpy.Reader.from_path(vcf_path)
    for record in reader:
        # Note: only take the first alternative
        outrow = [
            record.CHROM, record.POS, record.REF, record.ALT[0].serialize()
        ]
        filt = record.FILTER[0]  # PASS
        #print("INFO: ", record.INFO)
        wildtype = record.INFO.get("WT")
        het = record.INFO.get("HET")
        hom = record.INFO.get("HOM")
        nc = record.INFO.get("NC")

        #print("FORMAT: ", record.FORMAT)
        call = record.calls[0]
        #print(call)
        freq = call.data.get("FREQ")
        cov = call.data.get("SDP")
        outrow.extend([freq, cov])

        reads1 = call.data.get('RD')
        reads1_plus = call.data.get('RDF')
        reads1_minus = call.data.get('RDR')
        outrow.append("%-15s" % ("%s=%s:%s" % (reads1, reads1_plus, reads1_minus)))

        reads2 = call.data.get("AD")
        reads2_plus = call.data.get("ADF")
        reads2_minus = call.data.get("ADR")
        outrow.append("%-15s" % ("%s=%s:%s" % (reads2, reads2_plus, reads2_minus)))

        pval = call.data.get("PVAL")
        outrow.append("%-20s" % ("%s=%s" % (filt, pval)))
        outrow.extend([wildtype, het, hom, nc])

        print("\t".join(map(str, outrow)))


def process_varscan_file(infile):
    with open(args.infile) as f:
        for line in f:
            line = line.strip()
            a = line.split('\t')
            if a[0] == 'Chrom':
                continue
            # Take Chrom, Position, Ref, Var verbatim
            outrow = [a[0], a[1], a[2], a[3]]

            # split Cons:Cov:Reads1:Reads2:Freq:P-value
            b = a[4].split(':')

            # split StrandFilter:R1+:R1-:R2+:R2-:pval
            c = a[5].split(':')

            # add Freq, Cov
            outrow.extend([b[4], b[1]])

            if len(c) == 7:
                # This happens if there are numbered passes
                # In that case, make sure the Pass is copied together
                # with its number
                # example: Pass:1.0:84:75:2:0:2.8377E-1
                outrow.append("%-15s" % ("%s=%s:%s" % (b[2], c[2], c[3])))
                outrow.append("%-15s" % ("%s=%s:%s" % (b[3], c[4], c[5])))
                outrow.append("%-20s" % ("%s:%s=%s" % (c[0], c[1], c[6])))
                outrow.extend([a[6], a[7], a[8], a[9]])
            elif len(c) == 6:
                # normal case
                # Reads1=R1+:R1-
                outrow.append("%-15s" % ("%s=%s:%s" % (b[2], c[1], c[2])))
                # Reads2=R2+:R2-
                outrow.append("%-15s" % ("%s=%s:%s" % (b[3], c[3], c[4])))
                # StrandFilter=pval
                outrow.append("%-20s" % ("%s=%s" % (c[0], c[5])))
                # SamplesRef (WT), SamplesHet, SamplesHom, SamplesNC
                outrow.extend([a[6], a[7], a[8], a[9]])
            else:
                raise Exception("NO")
            print("\t".join(outrow))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('infile', help="input varscan/vcf file")

    args = parser.parse_args()
    if args.infile.endswith(".vcf"):
        process_vcf_file(args.infile)
    else:
        process_varscan_file(args.infile)
