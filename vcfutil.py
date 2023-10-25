#!/usr/bin/env python3

import os
import gzip

def is_file_empty(path):
    if path.endswith('.gz'):
        print("Found gzipped file\n")
        with gzip.open(path, 'rt') as infile:
            line = infile.readline()
            return len(line) == 0
    else:
        return os.stat(path).st_size == 0

def vcfopen(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _get_frequency_from_line(line, vcf_file_program):
    vector = line.split("\t")

    ## if vcf file is varscan file
    if vcf_file_program == "varscan_results":
        ## collect vcf field for frequency
        freq = vector[9].split(':')[6]

        #print('++++++++Freq: %s+++++++++' %(freq))
        chm = line.split('\t')[0]
        my_fields = vector[7].split(';')
        #adp = myFields.split("ADP=")[1]

        ## if vcf file is samtools file
    elif vcf_file_program == "samtools_results":
        #collect list of all fields from vcf file
        my_fields = vector[7].split(';')
        #print(my_fields)
        # grab DP4 field
        dp4 = next(filter(lambda x:'DP4' in x, my_fields))
        #print dp4
        dp4_fields = dp4.split("DP4=")[1]
        fr = float(dp4_fields.split(",")[0])
        rr = float(dp4_fields.split(",")[1])
        fa = float(dp4_fields.split(",")[2])
        ra = float(dp4_fields.split(",")[3])
        freq = (fa + ra) / (fr + rr + fa + ra)
        freq = round(freq * 100, 2)
        freq = str(freq) + "%"

        ## if vcf file is gatk file
    elif vcf_file_program == "gatk_results":
        ## collect vcf field
        my_fields = vector[9].split(':')
        #print myFields
        if my_fields[0] != ".":
            ad = float(my_fields[1].split(',')[1])
            dp = float(my_fields[2])
            freq = ad / dp
            freq = round(freq * 100, 2)
            freq = str(freq) + "%"
        else:
            freq = ''
    else:
        freq = ''
    return freq

def get_frequencies(snpeff_filtered_vcf, vcf_file_program):

    if is_file_empty(snpeff_filtered_vcf):
        print("FILE IS EMPTY !!!")
        return

    with vcfopen(snpeff_filtered_vcf) as infile:
        # skip all the comments
        prev_line = ''
        for line in infile:
            if not line.startswith("#"):
                break
            prev_line = line
        if not prev_line.startswith('#CHROM'):
            print("No chromosome information !!!, exiting combine_variants()")
            return

        frequencies = []
        # process the first non-comment line
        freq0 = _get_frequency_from_line(line, vcf_file_program)
        frequencies.append(freq0)
        for line in infile:
            frequencies.append(_get_frequency_from_line(line, vcf_file_program))

        print( 'Length of Frequencies= '+ str(len(frequencies)))
    return frequencies

####################### Combine variants from snpeff outputs ###############################
def combine_variants(snpeff_filtered_vcf, snpeff_final, combined_variants, runstats_outfile):
    """
    snpeff_filtered_vcf: input VCF file
    """
    print( "\033[34m Running Combine variants.. \033[0m")

    vcf_file_program = snpeff_filtered_vcf.split('/')[-2]
    print('snpeff_filtered_vcf is:%s' %snpeff_filtered_vcf)
    print('vcf_file_program:%s' %vcf_file_program)

    ## open vcf file for processing and skip comment lines
    frequencies = get_frequencies(snpeff_filtered_vcf, vcf_file_program)

    # filename for the output from converting vcf to oneliner with frequency and program added
    outfile_w_freq = combined_variants + '/' + snpeff_filtered_vcf.split('/')[-1].split('_snpeff_filtered')[0] + '_outfile_w_freq.txt'

    with open(outfile_w_freq, 'w') as outfile, open(snpeff_final, 'r') as final_in:
        final_in.readline()  # skip header

        myList = []
        index = 0
        for line in final_in:
            #if not line.strip():
            line = line.rstrip('\n')
            #print(line)
            # change alternative chromosome names
            chromosome = line.split('\t')[0]
            if chromosome == "NC_002937":
                print( "Found alternative chromosome name: %s" %chromosome)
                chm = "Chromosome"
            elif chromosome == "NC_005863":
                print( "Found alternative chromosome name: %s" %chromosome)
                chm = "pDV"
            else:
                chm = line.split('\t')[0]
            #print("Frequencies: " + frequencies[index] + "index: " + str(index))   
            pos = line.split('\t')[1]
            refs = list(line.split('\t')[2])[0] # grab only first character in reference sequence
            ref = line.split('\t')[2]
            alt = line.split('\t')[3]
            eff = line.split('\t')[9]
            imp = line.split('\t')[10]
            fnc = line.split('\t')[11]
            cdn = line.split('\t')[12]
            aac = line.split('\t')[13]
            loc = line.split('\t')[15]
            cod = line.split('\t')[16]
            dep = line.split('\t')[6]
            fre = frequencies[index]
            pro = vcf_file_program.split('_')[0]

            lineToWrite = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chm,pos,refs,ref,alt,eff,imp,fnc,cdn,aac,loc,cod,fre,pro,dep)
            outfile.write(lineToWrite)

            myList.append(lineToWrite)
            index += 1

    # write into final output file
    if runstats_outfile is not None:
        for element in myList:
            runstats_outfile.write('%s' %element)


#SNPEFF_FILTERED_VCF = 'vcf_testdata/cmyc_L1_samtools_snpeff_filtered.vcf.gz'
SNPEFF_FILTERED_VCF = 'vcf_testdata/snps.vcf.gz'
#SNPEFF_FILTERED_VCF = 'vcf_testdata/empty.vcf.gz'
SNPEFF_FINAL = 'vcf_testdata/cmyc_L1_samtools_snpeff_final.txt'
COMBINED_VARIANTS = 'vcf_testdata/combined_variants'

if __name__ == '__main__':
    combine_variants(SNPEFF_FILTERED_VCF, SNPEFF_FINAL, COMBINED_VARIANTS, None)
