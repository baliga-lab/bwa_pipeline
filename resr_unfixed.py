import vcfpy
import os, subprocess
import varscan as vs

GENOME_SIZE = 4411532.0
MIN_COVERAGE = 3


def get_result_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps_unfixed.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_result_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds_unfixed.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_filt_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_filt_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds_filtered.vcf' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_for_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.for' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_for_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.for' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_forup_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.forup' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_forup_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.forup' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mix_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.mix' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mix_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.mix' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixfor_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.mixfor' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixfor_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.mixfor' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixmark_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.mixmark' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixmark_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.mixmark' % vs.get_varscan_files_path(varscan_results, exp_name)

def get_mixmarkkept_snp_file(varscan_results, exp_name):
    return '%s_varscan_snps.mixmarkkept' % vs.get_varscan_files_path(varscan_results, exp_name)


def get_mixmarkkept_indel_file(varscan_results, exp_name):
    return '%s_varscan_inds.mixmarkkept' % vs.get_varscan_files_path(varscan_results, exp_name)


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
    cmd = ["java", "-jar", "/proj/omics4tb2/wwu/snpEff/snpEff.jar",
           "-ud", "0", "-classic", "-csvStats",
           "snpeff-indel-stat.txt", "-geneId", "-lof", "-v", "-formatEff",
           "-o", "gatk", "Mycobacterium_tuberculosis_h37rv",
           "SRR10040387-indel.final.vcf", ">",
           "SRR10040387-indel.snpeff-formatted.vcf"]
# snpsift
    #cat SRR10040387-indel.snpeff-formatted.vcf | java -jar /proj/omics4tb2/wwu/snpEff/SnpSift.jar filter -p "((FILTER = 'PASS') & (EFF[*].CODING != 'NON_CODING'))" > SRR10040387-indel.snpsift-filtered.vcf
    #cat SRR10040387-indel.snpsift-filtered.vcf | perl /proj/omics4tb2/wwu/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /proj/omics4tb2/wwu/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "GEN[*].AD" "GEN[*].RD" "(FILTER = 'PASS')" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].CODING" "EFF[*].RANK" "EFF[*].DISTANCE" > SRR10040387-indel.finalized.vcf

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
    run_snp_eff(varscan_results, exp_name, config)

