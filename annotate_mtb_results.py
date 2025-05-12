#!/usr/bin/env python3

import pandas
import vcfpy
import argparse
import os

MTB_GENOME_SIZE = 4411532

def read_databases(dbpath):
    # 1. read Tuberculist information
    header = ['GeneID', 'GeneName', 'Strand', 'Start', 'Stop', 'Description', 'Category']
    tuberculist_df = pandas.read_csv(os.path.join(dbpath, '2_Tuberculist_new_20150307'),
                                     sep='\t',
                                     header=None, names=header)
    gene_infos = {}
    gene_ids = []

    igr = {}
    igr_keys = []

    for index, row in tuberculist_df.iterrows():
        gene_ids.append(row['GeneID'])
        gene_infos[row['GeneID']] = {
            "name": row["GeneName"],
            "strand": row["Strand"],
            "start": row["Start"],
            "stop": row["Stop"],
            "description": row["Description"],
            "category": row["Category"]
        }
    # we have the gene_infos map, now build igr
    for index, row in tuberculist_df.iterrows():
        if index > 0:
            if row["Start"] > tuberculist_df.iloc[index - 1]["Stop"]:
                prev_gene = tuberculist_df.iloc[index - 1]["GeneID"]
                igr_key = "%s-%s" % (prev_gene, row["GeneID"])
                igr[igr_key] = {
                    "start": gene_infos[prev_gene]["stop"] + 1,
                    "stop": row['Start'] - 1
                }
                igr_keys.append(igr_key)

    # 2. read Codon -> AA mapping
    code = {}
    genetic_codes_df = pandas.read_csv(os.path.join(dbpath, '3_genetic_codes'),
                                       sep='\t',
                                       header=None, names=['Codon', 'Result'])
    for index, row in genetic_codes_df.iterrows():
        code[row['Codon']] = row['Result']

    # 3. Read genome
    genome = ""
    with open(os.path.join(dbpath, '4_mtbc_sequence.fasta')) as infile:
        for line in infile:
            if not line.startswith(">"):
                genome += line.strip()

    return gene_ids, gene_infos, code, genome, igr_keys, igr


def process_res_row(pos, alt, special_map, dbinfo):
    gene_ids, gene_infos, code, genome, igr_keys, igr = dbinfo
    special = special_map[pos]
    for gene_id in gene_ids:
        gene_info = gene_infos[gene_id]
        if is_in_gene(pos, gene_info):
            if gene_id.startswith('MTB'):
                print("%s\t%s\t%s\t-\t---\t---\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (pos, record.REF, alt,
                       gene_id,
                       gene_info['name'], gene_info['description'], gene_info['category'],
                       # varscan special
                       special["cons"],
                       special["cov"],
                       special["reads1"],
                       special["reads2"],
                       special["freq"]))
            else:
                if gene_info['strand'] == '+':
                    #print("NON-MTB GENE (+) %s\n" % gene_id);
                    length = gene_info['stop'] - gene_info['start'] + 1
                    seq = genome[gene_info['start'] - 1:(gene_info['start'] -  1) + length]
                    loci = pos - gene_info['start'] + 1
                    ct = int(loci / 3)
                    remain = loci % 3
                    if remain == 0:
                        codon = ct
                        wildtype = seq[loci - 3: (loci - 3) + 3]
                        mutation = wildtype[:2] + alt
                    elif remain == 1:
                        codon = ct + 1
                        wildtype = seq[loci - 1: (loci - 1) + 3]
                        mutation = alt + wildtype[1:3]
                    elif remain == 2:
                        codon = ct + 1
                        wildtype = seq[loci - 2: (loci - 2) + 3]
                        mutation = wildtype[0] + alt + wildtype[2]

                elif gene_info['strand'] == '-':
                    #print("NON-MTB GENE (-): %s\n" % gene_id);
                    length = gene_info['stop'] - gene_info['start'] + 1
                    seq = genome[gene_info['start'] - 1:(gene_info['start'] -  1) + length]
                    seq = seq[::-1]  # reverse the sequence
                    loci = gene_info['stop'] - pos + 1
                    ct = int(loci / 3)
                    remain = loci % 3
                    if remain == 0:
                        codon = ct
                        wildtype = seq[loci - 3: (loci - 3) + 3]
                        mutation = wildtype[0:2] + alt
                    elif remain == 1:
                        codon = ct + 1
                        wildtype = seq[loci - 1: (loci - 1) + 3]
                        mutation = alt + wildtype[1:3]
                    elif remain == 2:
                        codon = ct + 1
                        wildtype = seq[loci - 2: (loci - 2) + 3]
                        mutation = wildtype[0] + alt + wildtype[2]

                    # translate bases
                    wildtype = wildtype.translate(str.maketrans("ATGC", "TACG"))
                    mutation = mutation.translate(str.maketrans("ATGC", "TACG"))

                try:
                    if code[wildtype] == code[mutation]:
                        restype = "Synonymous"
                    else:
                        restype = "Nonsynonymous"
                except KeyError:
                    restype = "Nonsynonymous"

                try:
                    wt_aa = code[wildtype]
                except KeyError:
                    wt_aa = "(%s)" % wildtype
                try:
                    mut_aa = code[mutation]
                except KeyError:
                    mut_aa = "(%s)" % mutation

                code_info = '%s-%s-%s' % (restype, wt_aa, mut_aa)
                triplets = '%s_%s' % (wildtype, mutation)
                print("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (pos, record.REF, alt,
                       codon, code_info, triplets, gene_id, gene_info['name'],
                       gene_info['description'], gene_info['category'],
                       # varscan special
                       special["cons"],
                       special["cov"],
                       special["reads1"],
                       special["reads2"],
                       special["freq"]))

        else:
            #print("not a gene location")
            pass

    # II. Iterate over joined genes
    for j in igr_keys:
        j_info = igr[j]
        if pos >= j_info['start'] and pos <= j_info['stop']:
            gene1, gene2 = j.split('-')
            if gene1 == 'Rv3924c' and gene2 == 'Rv0001':
                left = pos - j_info['start']
                right = MTB_GENOME_SIZE - pos;
            else:
                left = pos - gene_infos[gene1]['stop']
                right= gene_infos[gene2]['start'] - pos
            gene_info1 = gene_infos[gene1]
            gene_info2 = gene_infos[gene2]
            loc_str = "%s%d-%d%s" % (gene_info1['strand'], left, right, gene_info2['strand'])
            name_str = "%s-%s" % (gene_info1['name'], gene_info2['name'])
            desc_str = "%s##%s" % (gene_info1['description'], gene_info2['description'])
            cat_str = "%s##%s" % (gene_info1['category'], gene_info2['category'])
            print("%d\t%s\t%s\t-\t---\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                  (pos, record.REF, alt,
                   loc_str, j, name_str, desc_str, cat_str,
                   # varscan special
                   special["cons"],
                   special["cov"],
                   special["reads1"],
                   special["reads2"],
                   special["freq"]))


def process_vcfreader(res_reader, special_map, dbpath):
    #gene_ids, gene_infos, code, genome, igr_keys, igr = read_databases(dbpath)
    dbinfo = read_databases(dbpath)
    header = ["VarscanPosition", "Ref", "Alt", "CodonPos", "Type_WTAA_MutAA",
              "WTCodon_MutCodon", "Gene ID", "Name",
              "Description", "Type", "Cons", "Cov", "Reads1", "Reads2", "Freq"]
    print('\t'.join(header))

    for record in res_reader:
        pos = record.POS
        alt = record.ALT[0].value
        process_res_row(pos, alt, special_map, dbinfo)


def process_df(df, special_map, dbpath):
    #gene_ids, gene_infos, code, genome, igr_keys, igr = read_databases(dbpath)
    dbinfo = read_databases(dbpath)
    header = ["VarscanPosition", "Ref", "Alt", "CodonPos", "Type_WTAA_MutAA",
              "WTCodon_MutCodon", "Gene ID", "Name",
              "Description", "Type", "Cons", "Cov", "Reads1", "Reads2", "Freq"]
    print('\t'.join(header))

    for index, row in df.iterrows():
        pos = row["POS"]
        alt = row["ALT"]
        process_res_row(pos, alt, special_map, dbinfo)


DESCRIPTION = """annotate_mtb_results - Annotate SNP results with Tuberculist information"""

def valueof(value):
    return "NA" if value is None else value


def is_in_gene(pos, gene_info):
    return pos >= gene_info['start'] and pos <= gene_info['stop']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('resfile', help="result VCF file (VCF result from resr unfixed)")
    parser.add_argument('dbpath', help="path for database files")
    parser.add_argument('vcffiles', nargs='+', help="original VCF file/s (SNP/INDEL)")
    parser.add_argument("--nonvcf", action="store_true", help="result file is not a VCF, but comes from fixed part")
    args = parser.parse_args()

    ##### BUILD SPECIAL MAP
    special_map = {}
    for vcffile in args.vcffiles:
        reader = vcfpy.Reader.from_path(vcffile)
        for record in reader:
            chrom = record.CHROM
            pos = record.POS
            #special = special_map[pos]
            call = record.calls[0]
            """
            print("GT PHASES: ", call.gt_bases, ", GT_PHASE CHAR: ",
                  call.gt_phase_char, ", GT ALLELES: ", call.gt_alleles,
                  ", GT TYPE: ", call.gt_type)
            """
            special_map[pos] = {
                "cov": valueof(record.INFO.get("ADP")),
                "cons": "%s" % valueof(call.data.get("GT")),  # genotype consensus
                "gt": valueof(call.gt_bases),
                "reads1": valueof(call.data.get("RD")),
                "reads2": valueof(call.data.get("AD")),
                "freq": valueof(call.data.get("FREQ")),
                "pval": valueof(call.data.get("PVAL"))
            }

    ################### END BUILDING SPECIAL MAP

    if args.nonvcf:
        res_df = pandas.read_csv(args.resfile, sep='\t', header=0)
        process_df(res_df, special_map, args.dbpath)
    else:
        res_reader = vcfpy.Reader.from_path(args.resfile)
        process_vcfreader(res_reader, special_map, args.dbpath)
