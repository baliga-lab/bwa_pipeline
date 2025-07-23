import vcfpy
import os, subprocess
import varscan as vs
from snpeff import SnpEff
import pandas
import traceback

GENOME_SIZE = 4411532.0
MIN_COVERAGE = 3


# The official file names for the RESR workflow are
# defined here
# Don't try to build them anywhere else, otherwise we will
# lose consistency

def get_result_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_result_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_snpeff_format_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_snps_snpeff_formatted.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_snpeff_format_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_inds_snpeff_formatted.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_snpsift_filtered_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_snps_snpsift_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_snpsift_filtered_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_inds_snpsift_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_finalized_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_snps_finalized.txt' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_finalized_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_inds_finalized.txt' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_annotated_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_snps_annotated.tsv' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_annotated_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_inds_annotated.tsv' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_filt_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_filt_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_for_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.for' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_for_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.for' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_forup_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.forup' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_forup_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.forup' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mix_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.mix' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mix_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.mix' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixfor_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.mixfor' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixfor_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.mixfor' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixmark_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.mixmark' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixmark_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.mixmark' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixmarkkept_snp_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_snps.mixmarkkept' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixmarkkept_indel_file(varscan_results, exp_name):
    return '%s_RESR_UNFIXED_varscan_inds.mixmarkkept' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_avg_sequencing_depth(cns_file):
    print("get_avg_sequencing_depth on:  ", cns_file)
    reader = vcfpy.Reader.from_path(cns_file)
    the_sum = 0
    num_records = 0
    for record in reader:
        adp = int(record.INFO.get('ADP'))
        if adp >= MIN_COVERAGE:
            the_sum += adp
            num_records += 1

    return (num_records / GENOME_SIZE, the_sum / num_records)


def drop_excluded(varscan_results, exp_name, config):
    print("drop_excluded() - filter out excluded loci")
    cns_file = vs.get_cns_file(varscan_results, exp_name)
    snp_file = vs.get_snps_file(varscan_results, exp_name)
    filt_snp_file = get_filt_snp_file(varscan_results, exp_name)
    indel_file = vs.get_indels_file(varscan_results, exp_name)
    filt_indel_file = get_filt_indel_file(varscan_results, exp_name)
    drop_excluded_script = os.path.join(config["run_dir"], "varscan_drop_excluded.py")

    cmd = [drop_excluded_script,
           os.path.join(config["resr_database_dir"], "Excluded_loci_mask.list"),
           indel_file, filt_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    cmd = [drop_excluded_script,
           os.path.join(config["resr_database_dir"], "Excluded_loci_mask.list"),
           snp_file, filt_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)


def format_transform(varscan_results, exp_name, config):
    print("format_transform() - transform filtered VCF to .for file")
    format_transform_script = os.path.join(config["run_dir"],
                                           "unfixed_format_trans.py")
    filt_snp_file = get_filt_snp_file(varscan_results, exp_name)
    filt_indel_file = get_filt_indel_file(varscan_results, exp_name)
    for_snp_file = get_for_snp_file(varscan_results, exp_name)
    for_indel_file = get_for_indel_file(varscan_results, exp_name)

    cmd = [format_transform_script, filt_indel_file, ">", for_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    cmd = [format_transform_script, filt_snp_file, ">", for_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)


def mix_pileup_merge(varscan_results, exp_name, config):
    print("mix_pileup_merge() - merge .for file with info from pileup")
    mix_pileup_merge_script = os.path.join(config["run_dir"],
                                           "mix_pileup_merge.py")
    for_snp_file = get_for_snp_file(varscan_results, exp_name)
    for_indel_file = get_for_indel_file(varscan_results, exp_name)
    pileup_file = vs.get_pileup_file(varscan_results, exp_name)
    forup_snp_file = get_forup_snp_file(varscan_results, exp_name)
    forup_indel_file = get_forup_indel_file(varscan_results, exp_name)

    cmd = [mix_pileup_merge_script, for_snp_file, pileup_file,
           ">", forup_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)
    cmd = [mix_pileup_merge_script, for_indel_file, pileup_file,
           ">", forup_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)


MIX_EXTRACT_MAX_FREQUENCY = 95
MIX_EXTRACT_REAL_READS_RATIO_THRESH = 0.8
MIX_EXTRACT_MIN_READS_THRESH = 4

def mix_extract(varscan_results, exp_name, config):
    print("mix_extract() - filter forup entries and transform to mix file")
    forup_snp_file = get_forup_snp_file(varscan_results, exp_name)
    mix_snp_file = get_mix_snp_file(varscan_results, exp_name)
    _mix_extract(forup_snp_file, mix_snp_file)

    forup_indel_file = get_forup_indel_file(varscan_results, exp_name)
    mix_indel_file = get_mix_indel_file(varscan_results, exp_name)
    _mix_extract(forup_indel_file, mix_indel_file)

def _mix_extract(forup_file, mix_file):
    with open(forup_file) as forup, open(mix_file, "w") as mix:
        for line in forup:
            line = line.strip()
            a = line.split("\t")
            mut_freq = float(a[4].replace("%", ""))
            if mut_freq <= MIX_EXTRACT_MAX_FREQUENCY:
                b = a[7].split("=")
                c = a[6].split("=")
                real = float(b[0]) + float(c[0])
                # 1. real-read ratio >= threshold
                # 2. min reads number >= threshold for both wildtype and mutant alleles
                if (real / float(a[5]) > MIX_EXTRACT_REAL_READS_RATIO_THRESH and
                    int(b[0]) >= MIX_EXTRACT_MIN_READS_THRESH and
                    int(c[0]) >= MIX_EXTRACT_MIN_READS_THRESH):
                    d = b[1].split(":")
                    e = c[1].split(":")
                    if int(d[0]) >= 1 and int(d[1]) >= 1 and int(e[0]) >= 1 and int(e[1]) >= 1:
                        mix.write("%s\n" % line)


def forup_format(varscan_results, exp_name, config):
    print("forup_format() - transform .mix format to .mixfor")
    mix_snp_file = get_mix_snp_file(varscan_results, exp_name)
    mix_indel_file = get_mix_indel_file(varscan_results, exp_name)
    mixfor_snp_file = get_mixfor_snp_file(varscan_results, exp_name)
    mixfor_indel_file = get_mixfor_indel_file(varscan_results, exp_name)

    forup_format_script = os.path.join(config["run_dir"],
                                       "forup_format.py")

    cmd = [forup_format_script, mix_snp_file, ">", mixfor_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    cmd = [forup_format_script, mix_indel_file, ">", mixfor_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

def info_mark(varscan_results, exp_name, config):
    mixfor_snp_file = get_mixfor_snp_file(varscan_results, exp_name)
    mixfor_indel_file = get_mixfor_indel_file(varscan_results, exp_name)
    mixmark_snp_file = get_mixmark_snp_file(varscan_results, exp_name)
    mixmark_indel_file = get_mixmark_indel_file(varscan_results, exp_name)

    info_mark_script = os.path.join(config["run_dir"], "info_mark.py")

    cmd = [info_mark_script, mixfor_snp_file, ">", mixmark_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)
    cmd = [info_mark_script, mixfor_indel_file, ">", mixmark_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

def redepin_filt(varscan_results, exp_name, config, ratio, depth):
    mixmark_snp_file = get_mixmark_snp_file(varscan_results, exp_name)
    mixmark_indel_file = get_mixmark_indel_file(varscan_results, exp_name)

    redepin_filt_script = os.path.join(config["run_dir"], "redepin_filt.py")

    cmd = [redepin_filt_script,
           os.path.join(config["resr_database_dir"], "Excluded_loci_mask.list"),
           str(ratio), str(depth),
           mixmark_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)
    cmd = [redepin_filt_script,
           os.path.join(config["resr_database_dir"], "Excluded_loci_mask.list"),
           str(ratio), str(depth),
           mixmark_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

def convert_to_vcf(varscan_results, exp_name, config):
    print("use .mixmarkkept to filter the original VCF file")
    filt_indel_file = get_filt_indel_file(varscan_results, exp_name)
    filt_snp_file = get_filt_snp_file(varscan_results, exp_name)

    mixmarkkept_snp_file = get_mixmarkkept_snp_file(varscan_results, exp_name)
    mixmarkkept_indel_file = get_mixmarkkept_indel_file(varscan_results, exp_name)
    result_snp_file = get_result_snp_file(varscan_results, exp_name)
    result_indel_file = get_result_indel_file(varscan_results, exp_name)

    filter_script = os.path.join(config["run_dir"], "filter_vcf_by_mixmark.py")
    cmd = [filter_script, filt_snp_file, mixmarkkept_snp_file,
           result_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)
    cmd = [filter_script, filt_indel_file, mixmarkkept_indel_file,
           result_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)


def run_snpeff(varscan_results, exp_name, config):
    print("Running snpEff")

    # Step 1: use snpEff to format the result
    print("1. using snpEff to format the result")
    result_indel_file = get_result_indel_file(varscan_results, exp_name)
    result_snp_file = get_result_snp_file(varscan_results, exp_name)
    snpeff_format_snp_file = get_snpeff_format_snp_file(varscan_results, exp_name)
    snpeff_format_indel_file = get_snpeff_format_indel_file(varscan_results, exp_name)
    snpsift_filtered_snp_file = get_snpsift_filtered_snp_file(varscan_results, exp_name)
    snpsift_filtered_indel_file = get_snpsift_filtered_indel_file(varscan_results, exp_name)

    cmd = [config["tools"]["snpeff"],
           "-ud", "0", "-classic", "-csvStats",
           os.path.join(varscan_results, "snpeff-indel-stat.txt"),
           "-geneId", "-lof", "-v", "-formatEff",
           "-o", "gatk", "Mycobacterium_tuberculosis_h37rv",
           result_indel_file, ">",
           # inds_snpeff_formatted.vcf
           snpeff_format_indel_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    cmd = [config["tools"]["snpeff"],
           "-ud", "0", "-classic", "-csvStats",
           os.path.join(varscan_results, "snpeff-snp-stat.txt"),
           "-geneId", "-lof", "-v", "-formatEff",
           "-o", "gatk", "Mycobacterium_tuberculosis_h37rv",
           result_snp_file, ">",
           # snps_snpeff_formatted.vcf
           snpeff_format_snp_file]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    # snpsift
    print("RUNNING SNPSIFT on '%s'" % snpeff_format_indel_file)
    cmd = [
        "cat",
        # inds_snpeff_formatted.vcf
        snpeff_format_indel_file,
        "|",
        config["tools"]["snpsift"],
        "filter", "-p",
        "\"((FILTER = 'PASS')\"",
        ">",
        # snps_snpsift_filtered.vcf
        snpsift_filtered_indel_file
    ]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    print("RUNNING SNPSIFT on '%s'" % snpeff_format_snp_file)
    cmd = [
        "cat",
        # snps_snpeff_formatted.vcf
        snpeff_format_snp_file,
        "|",
        config["tools"]["snpsift"],
        "filter", "-p",
        "\"((FILTER = 'PASS')\"",
        ">",
        # snps_snpsift_filtered.vcfsp
        snpsift_filtered_snp_file
    ]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    # SNPEFF

    print("RUNNING VOPL on '%s'" % snpsift_filtered_indel_file)
    print("RUNNING VOPL on '%s'" % snpsift_filtered_snp_file)
    snpeff = SnpEff(config['tools']['snpeff'], config['tools']['snpsift'], config['tools']['vceff_opl'])

    # vcfEffOnePerLine | snpSift


    # generate result in combined_variants format
    # NOT a VCF FILE !!! We annotate this differently
    finalized_indel_file = get_finalized_indel_file(varscan_results, exp_name)
    finalized_snp_file = get_finalized_snp_file(varscan_results, exp_name)
    print("EXTRACTING TO ", finalized_indel_file)
    print("EXTRACTING TO ", finalized_snp_file)

    snpeff.vceff_opl_extract_fields(snpsift_filtered_indel_file,
                                    finalized_indel_file,
                                    ["CHROM", "POS", "REF", "ALT", "AF",
                                     "\"(FILTER = \'PASS\')\"",
                                     "\"EFF[*].EFFECT\"",
                                     "\"EFF[*].IMPACT\"",
                                     "\"EFF[*].FUNCLASS\"",
                                     "\"EFF[*].CODON\"",
                                     "\"EFF[*].AA\"",
                                     "\"EFF[*].AA_LEN\"",
                                     "\"EFF[*].GENE\"",
                                     "\"EFF[*].CODING\"",
                                     "\"EFF[*].RANK\"",
                                     "\"EFF[*].DISTANCE\""]
                                    )

    # THE VCEFF_OPL_EXTRACT METHOD FAILS FOR THE allele frequency attribute, so we need to
    # extract the attributes by ourselves
    indel_attr_map = {}
    with open(snpsift_filtered_indel_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if line == "":
                continue
            row = line.split("\t")
            pos = row[1]
            af, vc, rc = extract_allele_freq_and_counts(row)
            indel_attr_map[pos] = (af, vc, rc)

    snpeff.vceff_opl_extract_fields(snpsift_filtered_snp_file,
                                    finalized_snp_file,
                                    ["CHROM", "POS", "REF", "ALT", "AF",
                                     "\"(FILTER = \'PASS\')\"",
                                     "\"EFF[*].EFFECT\"",
                                     "\"EFF[*].IMPACT\"",
                                     "\"EFF[*].FUNCLASS\"",
                                     "\"EFF[*].CODON\"",
                                     "\"EFF[*].AA\"",
                                     "\"EFF[*].AA_LEN\"",
                                     "\"EFF[*].GENE\"",
                                     "\"EFF[*].CODING\"",
                                     "\"EFF[*].RANK\"",
                                     "\"EFF[*].DISTANCE\""]
                                    )

    snp_attr_map = {}
    with open(snpsift_filtered_snp_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if line == "":
                continue
            row = line.split("\t")
            pos = row[1]
            af, vc, rc = extract_allele_freq_and_counts(row)
            snp_attr_map[pos] = (af, vc, rc)
    return indel_attr_map, snp_attr_map


def extract_allele_freq_and_counts(row):
    format_desc = row[8].split(":")
    data = row[9].split(":")
    freq_idx = format_desc.index("FREQ")
    vc_idx = format_desc.index("AD")
    rc_idx = format_desc.index("RD")
    return data[freq_idx].replace("%25", "%"), data[vc_idx], data[rc_idx]


def annotate_results(varscan_results, exp_name, config, indel_attr_map,
                     snp_attr_map):
    annot_script = os.path.join(config["run_dir"], "annotate_mtb_results.py")
    resr_snp_result = get_result_snp_file(varscan_results, exp_name)
    resr_indel_result = get_result_indel_file(varscan_results, exp_name)
    final_snp_result = get_finalized_snp_file(varscan_results, exp_name)
    final_indel_result = get_finalized_indel_file(varscan_results, exp_name)
    annotated_snp_file = get_annotated_snp_file(varscan_results, exp_name)
    annotated_indel_file = get_annotated_indel_file(varscan_results, exp_name)
    #print(indel_attr_map)
    #print(snp_attr_map)

    cmd = [
        annot_script,
        "--nonvcf",
        final_snp_result,
        config["resr_database_dir"],
        resr_snp_result,
        ">",
        annotated_snp_file
    ]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    cmd = [
        annot_script,
        "--nonvcf",
        final_indel_result,
        config["resr_database_dir"],
        resr_indel_result,
        ">",
        annotated_indel_file
    ]
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    # Last step: transform the VCF + annotated file to a final file
    # by combining the two VCF files into a large result file amended
    # by the Name column of the annotation file
    # 1. Make a pos2name map for everything
    pos2name = {}
    snp_df = pandas.read_csv(annotated_snp_file, sep="\t", header=0)
    indel_df = pandas.read_csv(annotated_indel_file, sep="\t", header=0)
    for index, row in snp_df.iterrows():
        pos = int(row["VarscanPosition"])
        name = row["Name"]
        pos2name[pos] = name

    for index, row in indel_df.iterrows():
        pos = int(row["VarscanPosition"])
        name = row["Name"]
        pos2name[pos] = name

    # TODO: FINAL SNP RESULT AND FINAL INDEL RESULT ARE NOT VCF
    # WE STILL HAVE TO CONCATENATE AND THEM INTO A RESULT
    snp_df = pandas.read_csv(final_snp_result, sep="\t", header=0)
    indel_df = pandas.read_csv(final_indel_result, sep="\t", header=0)

    with open(os.path.join(varscan_results, "RESR_FINAL_RESULTS.tsv"), "w") as out:
        header = [
            "CHROM", "POS", "1ST_REF_BASE", "ALT", "CHANGE", "EFFECT",
            "IMPACT", "CLASS",  "CODON", "AA_CHANGE",
            "GENE", "CODING",
            "Freq", "Reads1", "Reads2",
            "Name"
        ]
        out.write("\t".join(header) + '\n')

        for index, row in snp_df.iterrows():
            chrom = row["CHROM"]
            pos = int(row["POS"])
            name = pos2name.get(pos, "NA")
            alt = row["ALT"]
            out_row = [
                chrom, str(pos), row["REF"],
                alt[0],  # "alt"
                alt,  # change
                row["EFF[*].EFFECT"],
                row["EFF[*].IMPACT"],
                row["EFF[*].FUNCLASS"],
                row["EFF[*].CODON"],
                row["EFF[*].AA"],
                row["EFF[*].GENE"],
                row["EFF[*].CODING"],
                snp_attr_map[str(pos)][0],
                snp_attr_map[str(pos)][2],
                snp_attr_map[str(pos)][1],
                name
            ]
            out_row = list(map(str, out_row))
            out.write("\t".join(out_row) + "\n")

        for index, row in indel_df.iterrows():
            chrom = row["CHROM"]
            pos = int(row["POS"])
            name = pos2name.get(pos, "NA")
            alt = row["ALT"]  # TODO
            out_row = [
                chrom, pos,
                row["REF"],
                alt[0],  # alt
                alt,  # change
                row["EFF[*].EFFECT"],
                row["EFF[*].IMPACT"],
                row["EFF[*].FUNCLASS"],
                row["EFF[*].CODON"],
                row["EFF[*].AA"],
                row["EFF[*].GENE"],
                row["EFF[*].CODING"],
                indel_attr_map[str(pos)][0],
                indel_attr_map[str(pos)][2],
                indel_attr_map[str(pos)][1],
                name
            ]
            try:
                line = "\t".join(out_row) + "\n"
                out_row = list(map(str, out_row))
                out.write("\t".join(out_row) + "\n")
            except:
                traceback.print_exc()
                print("Error in position %s, skipping" % str(pos))


def run_resr_unfixed(varscan_results, exp_name, config):
    print("run_resr_unfixed()")
    cns_file = vs.get_cns_file(varscan_results, exp_name)
    all_ratio, avg_seq_depth = get_avg_sequencing_depth(cns_file)
    print("ALL_RATIO: ", all_ratio, " AVG SEQDEPTH: ", avg_seq_depth)

    # Run varscan_drop_excluded on snp and indel
    drop_excluded(varscan_results, exp_name, config)
    format_transform(varscan_results, exp_name, config)
    mix_pileup_merge(varscan_results, exp_name, config)
    mix_extract(varscan_results, exp_name, config)
    forup_format(varscan_results, exp_name, config)
    info_mark(varscan_results, exp_name, config)
    redepin_filt(varscan_results, exp_name, config, all_ratio, avg_seq_depth)
    convert_to_vcf(varscan_results, exp_name, config)
    indel_attr_map, snp_attr_map = run_snpeff(varscan_results, exp_name, config)
    annotate_results(varscan_results, exp_name, config, indel_attr_map, snp_attr_map)
