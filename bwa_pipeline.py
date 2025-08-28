#!/usr/bin/env python3
## Variant calling v2.0
import glob
import sys
import os
import re
from vcfutil import combine_variants
import vcfpy
import argparse
import json
import subprocess
import varscan
import pandas
import shutil

# use globalsearch's well tested and flexible
# pattern based FASTQ file discovery
from globalsearch.rnaseq.find_files import find_fastq_files
from samtools import SamTools
from bcftools import BcfTools
import varscan as vs
from gatk import GATK
from snpeff import SnpEff
import resr_unfixed


def run_bwa_index(genome_fasta, config):
    index_cmd = [config['tools']['bwa'], 'index', genome_fasta]
    compl_proc = subprocess.run(' '.join(index_cmd), shell=True, capture_output=False, check=True)

def create_genome_indexes(genome_fasta, config):
    # samtools fasta reference file indexing
    samtools = SamTools(config['tools']['samtools'])
    samtools.faidx(genome_fasta)

    gatk = GATK(config['tools']['gatk'])
    gatk.create_seq_dict(genome_fasta)

    run_bwa_index(genome_fasta, config)

############# Functions ##############

def create_dirs(samtools_results, gatk_results, varscan_results, data_trimmed_dir, fastqc_dir,
                alignment_results, combined_variants, tmp_dir):

    dirs = [samtools_results, gatk_results, varscan_results, data_trimmed_dir, fastqc_dir, alignment_results,
            combined_variants, tmp_dir]
    for dir in dirs:
        # create results folder
        print()
        print(dir)
        if not os.path.exists('%s' %(dir)):
            try:
                os.makedirs('%s' %(dir))
                print ('\033[31m %s directory does NOT exist. I am creating it. \033[0m' %(dir))
            except:
                print('WARNING: Something went wrong. I cannot create %s directory. Skipping.' %(dir))
        else:
            print ('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


####################### Trimgalore for quality and trimming ###############################

def run_trim_galore_single(fastqc_dir, data_trimmed_dir, fastq_file, config):
    print("\033[34m Running trim_galore (SINGLE) \033[0m")
    cmd = [config["tools"]["trim_galore"],
           "--fastqc_args", "\"--outdir %s\"" % fastqc_dir,  # FastQC argument
           "--output_dir", data_trimmed_dir,
           fastq_file]
    compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)


def run_trim_galore_paired(fastqc_dir, data_trimmed_dir, first_pair_file,
                           second_pair_file, config):
    print("\033[34m Running TrimGalore (PAIRED-END) \033[0m")
    cmd = [config["tools"]["trim_galore"],
           "--fastqc_args", "\"--outdir %s\"" % fastqc_dir,  # FastQC argument
           "--paired",
           "--output_dir", data_trimmed_dir,
           first_pair_file, second_pair_file]
    compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)


def trim_galore(first_pair_file, second_pair_file, folder_name, sample_id, file_ext,
               data_trimmed_dir, fastqc_dir, config):
    """run trim_galore with single/paired end selection and FastQC"""
    if second_pair_file is None:  # single
        run_trim_galore_single(fastqc_dir, data_trimmed_dir, first_pair_file, config)
    else:  # paired-end
        run_trim_galore_paired(fastqc_dir, data_trimmed_dir, first_pair_file,
                               second_pair_file, config)


def get_base_filename(alignment_results, sample_id):
    return os.path.join(alignment_results, sample_id)

def run_bwa_alignment(alignment_results, file_ext, first_file_name, second_file_name, lane,
                      folder_name, sample_id, RGId, RGSm, RGLb, RGPu, files_2_delete,
                      data_trimmed_dir, genome_fasta, config):
    """Run BWA for alignment"""
    print ("\033[34m Running BWA alignment... \033[0m")

    # define result files
    if second_file_name is None:  # single end reads
        second_pair_trimmed = None
        if file_ext == "gz":
            first_pair_trimmed = '%s/%s_trimmed.fq.gz' % (data_trimmed_dir, first_file_name)
        else:
            first_pair_trimmed = '%s/%s_trimmed.fq' % (data_trimmed_dir, first_file_name)
    else:
        if file_ext == "gz":
            first_pair_trimmed = '%s/%s_val_1.fq.gz'%(data_trimmed_dir,first_file_name)
            second_pair_trimmed = '%s/%s_val_2.fq.gz'%(data_trimmed_dir,second_file_name)
        else:
            first_pair_trimmed = '%s/%s_val_1.fq'%(data_trimmed_dir,first_file_name)
            second_pair_trimmed = '%s/%s_val_2.fq'%(data_trimmed_dir,second_file_name)

    print ('Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed, second_pair_trimmed))

    # modify read group information
    read_group = "'@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s\\tPU:%s'" %(RGId, RGSm, RGLb, RGPu)

    base_file_name = get_base_filename(alignment_results, sample_id)
    print('Base Filename: %s' %(base_file_name))

    # bwa run command
    if second_pair_trimmed is not None:
        cmd = "%s mem -t 16 -R %s %s %s %s > %s.sam" % (config['tools']['bwa'], read_group, genome_fasta,
                                                        first_pair_trimmed,
                                                        second_pair_trimmed,
                                                        base_file_name)
    else:
        # single read: run BWA with only one FASTQ file
        cmd = "%s mem -t 16 -R %s %s %s > %s.sam" % (config['tools']['bwa'], read_group, genome_fasta,
                                                     first_pair_trimmed,
                                                     base_file_name)

    print( "++++++ Run BWA Command", cmd)
    os.system(cmd)

    files_2_delete.append('%s.sam' % (base_file_name))
    return base_file_name


def run_samtools_fixmate_step(base_file_name, sample_id, files_2_delete, config):
    """samtools fixmate, sort and index to cleanup read pair info and flags"""
    print( "\033[34m Running SAMtools fixmate step... \033[0m")
    samtools = SamTools(config['tools']['samtools'])
    sorted_bam = "%s_sorted.bam" % base_file_name
    fixmate_bam = "%s_fixmate.bam" % base_file_name

    samtools.fixmate("%s.sam" % base_file_name, fixmate_bam)
    samtools.sort(fixmate_bam, sorted_bam, "%s/%s_temp" % (config["tmp_dir"], sample_id))
    samtools.index(sorted_bam)
    files_2_delete.extend([fixmate_bam, sorted_bam, '%s_sorted.bam.bai' % base_file_name])
    return sorted_bam


####################### GATK 1st Pass ###############################
def run_gatk(base_file_name,files_2_delete,exp_name,alignment_results):
    print("\033[34m Running GATK Realigner.. \033[0m")
    #run Target Interval Creater command java -Xmx128m -jar
    #cmd1 = '%s --java-options "-Xmx4g" -T RealignerTargetCreator -R %s -I %s_sorted.bam -o %s.intervals' % (GATK, genome_fasta, base_file_name, base_file_name)
    #run Indel Realigner command
    #cmd2 = '%s -T IndelRealigner -R %s -I %s_sorted.bam -targetIntervals %s.intervals -o %s_realigned.bam' %(GATK, genome_fasta, base_file_name, base_file_name, base_file_name)
    # index bam file
    #cmd3 = '%s index %s_sorted.bam' % (SAMTOOLS, base_file_name)
    # Detect covariates
    cmd4 = '%s BaseRecalibrator -R %s -I %s/%s_marked.bam --known-sites NA -O %s_recal.table' % (GATK, genome_fasta,alignment_results,exp_name,base_file_name)
    # Adjust quality scores
    cmd5 = '%s PrintReads -R %s -I %s/%s_marked.bam --BQSR %s_recal.table -O %s_recal.bam' % (GATK, genome_fasta, alignment_results,exp_name, base_file_name, base_file_name)

    #print("++++++ Command GATK Interval Creater: ", cmd1)
    #os.system(cmd1)
    #print( "++++++ Command GATK Realigner: ", cmd2)
    #os.system(cmd2)
    #print( "++++++ Command GATK index BAM: ", cmd3)
    #os.system(cmd3)
    #print("++++++ Command GATK BaseRecalibrator: ", cmd4)
    #os.system(cmd4)
    #print( "++++++ Command GATK PrintReads: ", cmd5)
    #os.system(cmd5)
    temp_files = ['%s.intervals'%base_file_name,'%s_realigned.bam'%base_file_name,'%s_realigned.bai'%base_file_name,'%s_recal.table'%base_file_name,'%s_recal.bam'%base_file_name,'%s_recal.bai'%base_file_name]
    for temp_file in temp_files:
        files_2_delete.append(temp_file)


def get_alignment_files_path(alignment_results, exp_name):
    return os.path.join(alignment_results, exp_name)


def mark_duplicates(alignment_results, exp_name, base_file_name, config):
    """Mark duplicates with Picard"""
    print("\033[34m Running Mark Duplicates.. \033[0m")
    print("Exp Name:", exp_name)
    # collect list of bwa aligned bam files
    aligned_bams = glob.glob('%s*_sorted.bam' % base_file_name)

    print('Output BAM File...')
    marked_bam_name = '%s/%s_marked.bam' % (alignment_results, exp_name)
    metrics_file = '%s/%s.metrics' % (alignment_results, exp_name)
    alignment_files_path = get_alignment_files_path(alignment_results, exp_name)

    picard_cmd = config['tools']['picard']
    # Mark Duplicates, use new format for picard command line
    markdup_cmd = [picard_cmd, "MarkDuplicates"]
    # create command line parameter for each file
    for aligned_bam in aligned_bams:
        markdup_cmd += ["-INPUT", aligned_bam]
    markdup_cmd += ["-VALIDATION_STRINGENCY", "LENIENT",
                    "-M", metrics_file,
                    "-O", marked_bam_name,
                    # include the sort order so Picard knows what to do
                    "-ASO", "coordinate"]

    print( "++++++ Mark Duplicated Command: '%s'... " % ' '.join(markdup_cmd))
    compl_proc = subprocess.run(' '.join(markdup_cmd), shell=True, capture_output=False, check=True)

    # index BAM file with Samtools
    samtools = SamTools(config['tools']['samtools'])
    samtools.index(marked_bam_name)


####################### Samtools Variant Calling ###############################
def bcftools_variants(samtools_results, alignment_files_path, exp_name, folder_name, config):
    print("\033[34m Running bcftools Variant Calling.. \033[0m")

    # create samtools results specific results directory
    # WW: This is a duplicate TODO REMOVE
    samtools_files_path = '%s/%s' % (samtools_results, exp_name)

    # Produce BCF file with all locations in the genome
    bcftools = BcfTools(config['tools']['bcftools'], config['tools']['tabix'],
                        config['tools']['bgzip'])

    # Prepare vcf file for querying with tabix. WW: This seems to be incomplete
    #cmd2 = '%s -p vcf %s_samtools.vcf.gz' % (TABIX, samtools_files_path)

    vcf_path = '%s_samtools.vcf' % samtools_files_path
    vcf_path2 = '%s_samtools2.vcf' % samtools_files_path
    vcf_path_final = '%s_samtools_final.vcf' % samtools_files_path
    bcftools.variant_calling_mpileup(genome_fasta, '%s_marked.bam' % alignment_files_path,
                                     vcf_path)
    bcftools.filter_variants(vcf_path, vcf_path2)
    bcftools.view(vcf_path2, vcf_path_final)


####################### Varscan Variant Calling ###############################
def varscan_variants(alignment_files_path, varscan_results, exp_name, folder_name, files_2_delete, config): # with varscan
    print("\033[34m Running Varscan.. \033[0m")

    # create varscan results specific results directory
    #varscan_files_path = '%s/%s' % (varscan_results, exp_name)
    pileup_file = vs.get_pileup_file(varscan_results, exp_name)

    # samtools mpileup
    samtools = SamTools(config['tools']['samtools'])
    samtools.mpileup(genome_fasta, "%s_marked.bam" % alignment_files_path,
                     pileup_file)

    varscan = vs.VarScan(config)
    varscan.mpileup2snp(varscan_results, exp_name)
    varscan.mpileup2indel(varscan_results, exp_name)
    varscan.mpileup2cns(varscan_results, exp_name)

    files_2_delete.append(pileup_file)


def gatk_variants(alignment_files_path, gatk_results, exp_name, folder_name, config):
    """GATK Variant Calling with HaploTypeCaller"""

    print("\033[34m Running GATK Haplotype Variant Caller.. \033[0m")
    # create varscan results specific results directory
    gatk_files_path = '%s/%s' % (gatk_results, exp_name)

    gatk = GATK(config['tools']['gatk'])
    marked_bam = '%s_marked.bam' % alignment_files_path
    raw_vcf = '%s_gatk_raw.vcf' % gatk_files_path
    gatk.haplotype_caller(genome_fasta, marked_bam, raw_vcf)

    # Select SNP variants
    snp_vcf = '%s_gatk_snps.vcf' % gatk_files_path
    filtered_snp_vcf = '%s_gatk_snps_filtered.vcf' % gatk_files_path
    gatk.select_variants(genome_fasta, raw_vcf, snp_vcf, "SNP")
    gatk.variant_filtration(genome_fasta, snp_vcf, filtered_snp_vcf,
                            "my_snp_filter",
                            "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0")

    # Select indel variants
    indel_vcf = "%s_gatk_inds.vcf" % gatk_files_path
    filtered_indel_vcf = "%s_gatk_inds_filtered.vcf" % gatk_files_path
    gatk.select_variants(genome_fasta, raw_vcf, indel_vcf, "INDEL")
    gatk.variant_filtration(genome_fasta, indel_vcf, filtered_indel_vcf,
                            "my_indel_filter",
                            "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0")

    # Merge vcf files
    cmd6 = "%s MergeVcfs I=%s_gatk_snps_filtered.vcf I= %s_gatk_inds_filtered.vcf O=%s_gatk_final.vcf" % (config['tools']['picard'], gatk_files_path, gatk_files_path, gatk_files_path)
    print("++++++ Merging vcf files: ", cmd6)
    os.system(cmd6)


####################### Run SNPEff annotations ###############################
def run_snpeff(samtools_files_path, varscan_files_path, gatk_files_path, combined_variants, exp_name,
               snpeff_db, config):
    print("\033[34m Running SNPEff Annotations.. \033[0m")

    vcf_files = ['%s_samtools_final.vcf' % samtools_files_path,
                 '%s_varscan_inds_final.vcf' % varscan_files_path,
                 '%s_varscan_snps_final.vcf' % varscan_files_path,
                 '%s_gatk_final.vcf' % gatk_files_path]

    # VCF file list to process
    print('input vcf_files:\n')
    for vcf_file in vcf_files:
        print(vcf_file + '\n')

    #varscan_files = ['%s_varscan_inds_final.vcf'%(varscan_files_path), '%s_varscan_snps_final.vcf'%(varscan_files_path)]
    #print('varscan_files:%s' %(varscan_files))

    # create output file for combined variants output
    combined_variants_output = '%s/%s_combined_variants.txt' %(combined_variants,exp_name)
    print("combined_variants_output: '%s'" %(combined_variants_output))

    # open the final output file for writing before combining
    run_snpeff.t = open(combined_variants_output, 'w')

    for vcf_file in vcf_files:
        print()
        print("Processing VCF file:" + vcf_file + '\n')

        ## WW: is this the correct comparison ?? Note the ".gz", but
        ## bcftools returns an uncompressed file
        if vcf_file == '%s_samtools_final.vcf.gz' % (samtools_files_path):
            # reformat file for gzip compression
            #plain_vcf_name = re.split('.gz', vcf_file)[0]
            #os.system('mv %s %s' %(vcf_file, plain_vcf_name))
            #os.system('%s view -Oz -o %s %s' % (BCFTOOLS, vcf_file, plain_vcf_name))
            #os.system('%s index %s' % (BCFTOOLS, vcf_file))

            # Rename chromosome names to match snpeff database
            vcf_file_renamed = re.split('final.', vcf_file)[0] + 'renamed.' + re.split('final.', vcf_file)[1]
            print("vcf_file_renamed: " + vcf_file_renamed)

            # gzip and index vcf file with tabix
            bgzip_cmd1 = '%s -c %s > %s' % (config['tools']['bgzip'], vcf_file, vcf_file)
            os.system(bgzip_cmd1)

            tabix_cmd1 = '%s -p vcf %s' % (config['tools']['tabix'], vcf_file)
            os.system(tabix_cmd1)

            bcft_cmd1 = '%s annotate --rename-chrs %s/chrom_names.txt %s > %s' % (config['tools']['bcftools'],
                                                                                  genome_dir,vcf_file,vcf_file_renamed)
            print("Renaming chromosomes for snpEFF with bcftools:" + bcft_cmd1)
            os.system(bcft_cmd1)

        else:
            # Rename chromosome names to match snpeff database
            vcf_file_renamed = re.split('final.', vcf_file)[0] + 'renamed.' + re.split('final.', vcf_file)[1]
            print("vcf_file_renamed: " + vcf_file_renamed)

            vcf_file_bgzip = re.split('final.', vcf_file)[0] + 'renamed.' + re.split('final.', vcf_file)[1] + ".bgz"

            # gzip and index vcf file with tabix
            bgzip_cmd1 = '%s -c %s > %s' % (config['tools']['bgzip'], vcf_file, vcf_file_bgzip)
            os.system(bgzip_cmd1)

            tabix_cmd1 = '%s -p vcf %s' % (config['tools']['tabix'], vcf_file_bgzip)
            os.system(tabix_cmd1)

            bcft_cmd1 = '%s annotate --rename-chrs %s/chrom_names.txt %s > %s' % (config['tools']['bcftools'],
                                                                                  genome_dir,
                                                                                  vcf_file_bgzip, vcf_file_renamed)
            print("Renaming chromosomes for snpEFF with bcftools:" + bcft_cmd1)
            os.system(bcft_cmd1)

        # create file names for output
        snpeff_vcf = re.split('final.', vcf_file)[0] + 'snpeff.' + re.split('final.', vcf_file)[1]
        print('snpeff_vcf:' + snpeff_vcf)

        snpeff_filtered_vcf = re.split('renamed.', vcf_file_renamed)[0] + 'snpeff_filtered.' + re.split('renamed.', vcf_file_renamed)[1]
        print('snpeff_filtered_vcf:' + snpeff_filtered_vcf)
        snpeff_stats = re.split('final.', vcf_file)[0] + 'snpeff_stats.txt'
        snpeff_final = re.split('final.', vcf_file)[0] + 'snpeff_final.txt'

        snpeff = SnpEff(config['tools']['snpeff'], config['tools']['snpsift'], config['tools']['vceff_opl'])

        # Run the snpEff commands
        snpeff.format_eff(snpeff_stats, snpeff_db, vcf_file_renamed, snpeff_vcf)
        snpeff.snpsift_filter(snpeff_vcf, snpeff_filtered_vcf,
                              "((FILTER = \'PASS\') & (EFF[*].CODING != \'NON_CODING\'))")
        snpeff.vceff_opl_extract_fields(snpeff_filtered_vcf, snpeff_final,
                                        ["CHROM", "POS", "REF", "ALT", "AF", "AC", "DP", "MQ",
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
                                         "\"EFF[*].DISTANCE\""])

        # Run function to combine vcf files into a single file from 3 callers
        combine_variants(snpeff_filtered_vcf, snpeff_final, combined_variants, run_snpeff.t)
    run_snpeff.t.close()

####################### Collate variants from three programs ###############################


def read_ppe_ins_loci(config):
    result = set()
    with open(os.path.join(config["resr_database_dir"], "PPE_INS_loci.list")) as f:
        for line in f:
            result.add(line.strip())
    return result

def salomon(aggregatedList, config):
    # this function takes a list of lists of variants and find the consensus list
    if config['organism'] == 'mtb':
        is_mtb = True
        ppe_ins_loci = read_ppe_ins_loci(config)
    else:
        is_mtb = False
        ppe_ins_loci = None

    # f.1. finding the unique positions
    uniqueLocations = []
    for variants in aggregatedList:
        for variant in variants:
            uniqueLocation = variant[:3]
            # special case for MTB: Exclude ppe_ins_loci
            if is_mtb:
                if uniqueLocation[1] in ppe_ins_loci:
                    continue

            if uniqueLocation not in uniqueLocations:
                uniqueLocations.append(uniqueLocation)

    # f.2. building the full consensus list
    consensus_list=[]
    for uniqueLocation in uniqueLocations:
        callers=[]
        body=[]
        freqs = []
        dps = []
        freq=''
        dp =''
        freqFloats = []

        for variants in aggregatedList:
            for variant in variants:
                if uniqueLocation == variant[:3]:
                    body = variant[:-2]

                    #if variant[-2] == 'varscan':
                    #    print("VARSCAN, BODY START: ", body)

                    callers.append(variant[-2])
                    if variant[-2] == 'varscan':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        #print variant[:-2]
                        dp = "NA"
                        dps.append(dp)

                        if freq == '':
                            print( 'WARNING varscan did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'samtools':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        dp =variant[-1]
                        dps.append(dp)
                        if freq == '':
                            print( 'WARNING samtools did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'gatk':
                        dp = variant[-1]
                        dps.append(dp)
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)

        # incorporating what we found
        callersString=':'.join(callers)
        freqsString=':'.join(freqs)
        dpsString=':'.join(dps)
        freqFloatString = max(freqFloats)

        # Get the samtools frequency as final frequency if available otherwise pick varscan freq.
        if 'samtools' in callers and len(callers) > 1:
            final_freq = freqs[callers.index('samtools')]
        elif 'varscan' in callers and 'samtools' not in callers:
            final_freq = freqs[callers.index('varscan')]
        else:
            final_freq = freqs[0]

        body.append(callersString)
        body.append(freqsString)
        body.append(dpsString)
        body.append(freqFloatString)
        body.append(final_freq)
        consensus_list.append(body)

    return consensus_list


####################### Variant Retriever ###############################
def variantRetriever(combined_output_file):

    # this function retrieves the variants for each caller
    varscan_list = []
    gatk_list = []
    samtools_list = []
    print('combined_output')
    print(os.stat(combined_output_file).st_size)

    # Start reading each line and append to appropriate program list
    with open(combined_output_file) as f:
        for line in f:
            vector = line.split('\t')
            vector[-1]=vector[-1].replace('\n','')
            program=vector[-2]
            #print('program is:%s' %(program))
            if program == 'varscan':
                varscan_list.append(vector)
            if program == 'gatk':
                gatk_list.append(vector)
            if program == 'samtools':
                samtools_list.append(vector)
        f.close()

    return varscan_list,gatk_list,samtools_list


####################### Collate variants ###############################
def collate_variants(exp_name, combined_output_file, merged_variants_file,
                     combined_variants, varscan_results,
                     gatk_files_path, samtools_files_path,
                     config):
    print( "\033[34m Running Collate variants.. (Combined output file: '%s', merged variants file: '%s') \033[0m" % (combined_output_file, merged_variants_file))

    # 2. recovering list of variants
    varscan_list,gatk_list,samtools_list = variantRetriever(combined_output_file)
    print( 'detected variants',len(varscan_list),len(gatk_list),len(samtools_list))

    # 3. finding consensus list of variants
    print( 'merging...')
    consensus_list = salomon([varscan_list, gatk_list, samtools_list], config)
    print( 'final set',len(consensus_list))

    # 4. writing a consensus list
    print( 'writing file...')

    with open(merged_variants_file, 'w') as g:
        g.write('CHROM\tPOS\t1ST_REF_BASE\tALT\tCHANGE\tEFFECT\tIMPACT\tCLASS\tCODON\tAA_CHANGE\tGENE\tCODING\tSAMTOOLS_FREQ\tVARIANT_CALLERS\tVARIANT_FREQS\tVARIANT_READS\tFREQ\tF1\n')
        for element in consensus_list:
            line2write='\t'.join(element)
            line2write=line2write+'\n'
            g.write(line2write)

    if config['organism'] == 'mtb':
        annotate_combined_results(exp_name, merged_variants_file, combined_variants,
                                  varscan_results,
                                  gatk_files_path, samtools_files_path,
                                  config)

def annotate_combined_results(exp_name, merged_variants_file, combined_variants,
                              varscan_results,
                              gatk_files_path, samtools_files_path,
                              config):
    """Tuberculist annotation and create a new result file based on the TSV file with a 'Name' column which
    has the intergenic relationship"""

    # Step 1. create annotation script based on the merged_variants_file and name
    # in the variable annotated_result
    annot_script = os.path.join(config["run_dir"], "annotate_mtb_results.py")
    all_snps_vcf = varscan.get_snps_file(varscan_results, exp_name)
    all_indels_vcf = varscan.get_indels_file(varscan_results, exp_name)
    annotated_result = os.path.join(combined_variants, "COMBINED_ANNOTATED.tsv")
    final_gatk_vcf = "%s_gatk_final.vcf" % gatk_files_path
    final_samtools_vcf = '%s_samtools_final.vcf' % samtools_files_path

    cmd = [
        annot_script,
        "--nonvcf",
        merged_variants_file,
        config["resr_database_dir"],
        all_snps_vcf, all_indels_vcf,
        final_gatk_vcf, final_samtools_vcf,
        ">",
        annotated_result
    ]
    print(' '.join(cmd))
    proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False,
                          check=True)

    # Step 2. use the "Name" column of the annotated_result and append it to the end
    # make a map from position to Name
    pos2name = {}
    df1 = pandas.read_csv(annotated_result, sep='\t', header=0)
    for index, row in df1.iterrows():
        pos2name[int(row["VarscanPosition"])] = row["Name"]

    # of the merged_variants file, aligned by the position
    tmp_merged = os.path.join(combined_variants, "tmp_annotated_merged.tsv")
    with open(merged_variants_file) as infile, open(tmp_merged, "w") as outfile:
        header = list(infile.readline().strip().split("\t"))
        header.append("Name")
        outfile.write("\t".join(header) + "\n")
        for line in infile:
            comps = line.strip().split("\t")
            pos = int(comps[1])
            name = pos2name.get(pos, "NA")
            comps.append(name)
            outfile.write("\t".join(comps) + "\n")

    # Step 3. move the temp file to the final file with the amendment
    shutil.move(tmp_merged, merged_variants_file)


def delete_temp_files(files_2_delete):
    print
    print( "\033[34m Deleting Temporry files.. \033[0m")
    for file in files_2_delete:
        cmd = 'rm %s' %(file)
        print(cmd)
        #os.system(cmd)


def run_tbprofiler(sorted_bam_file, sample_id, result_dir, config):
    """tb-profiler profile --bam ERR2652935.sralite.1_1_sorted.bam --dir tbprofileout --prefix ERR2652935 --csv"""
    # activate TBProfiler environment
    profile_cmd = [
        # we have to run tb profiler in a separate conda environment !!!
        "conda", "run", "-n", config["tbprofiler_env"], "--live-stream",
        config["tools"]["tbprofiler"], "profile", "--bam", sorted_bam_file,
        "--dir", result_dir,
        "--prefix", sample_id, "--csv"]
    print("Running TBprofiler command: '%s'" % " ".join(profile_cmd))
    compl_proc = subprocess.run(" ".join(profile_cmd),
                                shell=True, capture_output=False, check=True)


def run_pipeline(organism, data_folder, resultdir, snpeff_db, genome_fasta, config):
    print("run_pipeline()")
    # retrieve list of pairs. paired end data has both elements set,
    # single read has the second element of the pair set to None
    fastq_files = find_fastq_files(data_folder, config['fastq_patterns'])
    print(fastq_files)

    folder_count = 1
    files_2_delete = [] # create a list of files to delete later and keep adding to the list

    # Get the folder name
    folder_name = data_folder.split('/')[-1]
    print( '\033[33mProcessing Folder: %s \033[0m' % folder_name)

    # get the list of first file names in paired end sequences
    first_pair_files = [f[0] for f in fastq_files]
    second_pair_files = [f[1] for f in fastq_files if f[1] is not None]

    print('There are %s fastq files in this directory.' % (len(first_pair_files) * 2))
    print("First pair files:")
    for i in first_pair_files:
        print(i)

    print("Second pair files:")
    for i in second_pair_files:
        print(i)

    # Program specific results directories
    folder_results_dir = "%s/%s/%s" % (resultdir, organism, folder_name)
    samtools_results = os.path.join(folder_results_dir, "samtools_results")
    gatk_results = os.path.join(folder_results_dir, "gatk_results")
    varscan_results = os.path.join(folder_results_dir, "varscan_results")
    alignment_results = os.path.join(folder_results_dir, "alignment_results")
    combined_variants = os.path.join(folder_results_dir, "combined_variants")
    tbprofiler_results = os.path.join(folder_results_dir, "tbprofiler")
    tmp_dir = os.path.join(config["tmp_dir"])

    # final results files
    combined_output_file = os.path.join(combined_variants, '%s_combined_variants.txt' % folder_name)
    print("combined_output_file: '%s'" % combined_output_file)
    merged_variants_file = os.path.join(combined_variants, '%s_merged_variants_final.txt' % folder_name)
    print("combined_output_file: '%s'" % combined_output_file)

    exp_name = folder_name
    alignment_files_path = get_alignment_files_path(alignment_results, exp_name)
    samtools_files_path = '%s/%s' % (samtools_results, exp_name)
    varscan_files_path = vs.get_varscan_files_path(varscan_results, exp_name)
    gatk_files_path = '%s/%s' % (gatk_results, exp_name)

    # 00. Get directories
    create_dirs(samtools_results, gatk_results, varscan_results,
                data_trimmed_dir, fastqc_dir, alignment_results,
                combined_variants, tmp_dir)

    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file, second_pair_file in fastq_files:
        first_file_name_full = os.path.basename(first_pair_file)
        file_ext = first_pair_file.split('.')[-1]
        second_file_name_full = os.path.basename(second_pair_file) if second_pair_file is not None else None

        print('\033[32m Processing Pair: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))
        first_file_name = re.split('.fastq|.fastq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fastq|.fastq.gz', second_file_name_full)[0] if second_file_name_full is not None else None
        print("first_file_name: '%s', second_file_name: '%s'" % (first_file_name, second_file_name))

        # Collect Sample attributes. The original extracts the a for A_B_C_L001_R1_001
        # and takes "C" for the lane
        try:
            lane = first_file_name.split("/")[-1].split("_")[2]
        except:
            lane = "SNA"  # not available in file name
        print("Lane: %s" %(lane))
        sample_id = re.split('.fastq|.fastq.gz', first_file_name)[0]
        print("sample_id: %s"  %(sample_id))

        # create readgroup info
        RGId = sample_id
        RGSm = exp_name
        RGLb = exp_name + "_l"
        RGPu = lane
        print( "RG: ID: %s SM: %s LB: %s PU: %s" %(RGId, RGSm, RGLb, RGPu))

        # Get the base file name here, because on reruns, checkpointing will leave
        # base_file_name undefined if not retrieved before the check
        base_file_name = get_base_filename(alignment_results, sample_id)
        # 01. Run trim_galore + FastQC
        trim_galore(first_pair_file, second_pair_file, folder_name, sample_id, file_ext,
                    data_trimmed_dir, fastqc_dir, config)

        # 02. Run bwa alignment to produce SAM file.
        run_bwa_alignment(alignment_results, file_ext, first_file_name,
                          second_file_name, lane,folder_name, sample_id,
                          RGId, RGSm, RGLb, RGPu, files_2_delete,
                          data_trimmed_dir, genome_fasta, config)

        # 03. Run samtools fixmate
        sorted_bam_file = run_samtools_fixmate_step(base_file_name, sample_id,files_2_delete, config)

        # 03b. Run TB Profiler on the sorted bam file
        if organism == 'mtb':
            try:
                tbprof_tool = config["tools"]["tbprofiler"]
                ## TODO: comment me in skipping for speed
                run_tbprofiler(sorted_bam_file, sample_id, tbprofiler_results, config)
            except:
                print("WARNING: can't find TBprofiler setting, skipping")

        file_count += 1

    # 05. Run Mark duplicates
    # WW: base_file_name is dependent on the loop to be run, does that make any sense ?
    mark_duplicates(alignment_results, exp_name, base_file_name, config)

    # 04. Run GATK 1st PASS
    run_gatk(base_file_name,files_2_delete,exp_name,alignment_results)

    # 06. Run Samtools variant calling
    bcftools_variants(samtools_results, alignment_files_path, exp_name, folder_name, config)
    # 07. Run varscan variant calling
    # TODO: We should look into this, because we call varscan 3 times,
    # but we actually might only need to call the mpileup2cns function
    # and extract the indels and snps results from the cns file
    varscan_variants(alignment_files_path, varscan_results, exp_name,
                     folder_name, files_2_delete, config)

    # 07b. Run resr unfixed pipeline
    resr_unfixed.run_resr_unfixed(varscan_results, exp_name, config)

    # 08. Run GATK variant Calling
    gatk_variants(alignment_files_path, gatk_results, exp_name, folder_name, config)

    # 09. Run SNPEff annotations
    vcf_file = run_snpeff(samtools_files_path, varscan_files_path, gatk_files_path, combined_variants, exp_name,
                          snpeff_db, config)
    # 10. Collate variants into single file from 3 callers and unify them
    collate_variants(exp_name, combined_output_file, merged_variants_file,
                     combined_variants, varscan_results,
                     gatk_files_path, samtools_files_path,
                     config)

    # 11. Delete temporary files_2_delete
    #delete_temp_files(files_2_delete)

    folder_count += 1


DESCRIPTION = "bwa_pipeline.py - BWA pipeline V2.0"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('config', help="configuration json file")
    # This is likely something like the sample number
    parser.add_argument('input_folder', help="input folder without the path")
    parser.add_argument('--index_genome',
                        help="run BWA indexing", action="store_true")

    # run dir needs to contain "reference" and will have a results directory
    args = parser.parse_args()
    with open(args.config) as infile:
        config = json.load(infile)

    ############# Data and Results directories ##############
    if not os.path.exists(config["result_dir"]):
        os.makedirs(config["result_dir"])

    # First set the top level analysis directory
    data_folder = '%s/%s' % (config["data_dir"], args.input_folder)

    # trimmed dir is in results because you can't assume you can write
    # to data folder
    organism = config["organism"]
    data_trimmed_dir = "%s/trimmed" % config["result_dir"]
    fastqc_dir = "%s/%s/%s/fastqc_results" % (config["result_dir"],
                                              organism,
                                              args.input_folder)

    known_sites = '%s-variants-compiled_sorted.vcf' % organism
    genome_dir = config["genome_dir"]
    print("genome dir: '%s'" % genome_dir)
    if not os.path.exists(genome_dir):
        exit("Genome directory: '%s' does not exist !!!" % genome_dir)
    try:
        org_config = config['organisms'][organism]
    except:
        exit("Organism '%s' not recognized !!!")

    ######### Annotation databases ############################
    # snpEff databases
    genome_gff = glob.glob(os.path.join(genome_dir, org_config["genome_gff"]))
    genome_fasta_path = os.path.join(genome_dir, org_config["genome_fasta"])
    try:
        print("GENOME FASTA: '%s'" % genome_fasta_path)
        genome_fasta = glob.glob(genome_fasta_path)[0]
    except:
        exit("Genome FASTA file '%s' not found, please check path" % genome_fasta_path)

    if args.index_genome:
        # create genome indexes. This should be optional
        create_genome_indexes(genome_fasta, config)
    else:
        snpeff_db = config["organisms"][organism]['snpeff_db']
        run_pipeline(organism, data_folder, config["result_dir"],
                     snpeff_db, genome_fasta, config)
