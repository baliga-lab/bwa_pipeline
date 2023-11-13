"""gatk.py - Interface to the gatk command line tool

This module provides a Python interface to the gatk command
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class GATK:
    """This interface was created based on GATK 4.3.0.0"""

    def __init__(self, cmd_path):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.cmd = cmd_path

    def create_seq_dict(self, genome_fasta):
        createdict_cmd = [self.cmd, "CreateSequenceDictionary", "-R", genome_fasta]
        # make the dict file name
        genome_name = genome_fasta.replace(".fna", "").replace(".fasta", "")
        #genome_name = genome_fasta.split(".fna")[0]
        #genome_name = genome_fasta.split(".fna")[0]
        if not os.path.exists('%s.dict' % genome_name):
            print ('\033[31m %s.dict  DOES NOT exist, creating. \033[0m' % genome_fasta)
            compl_proc = subprocess.run(' '.join(createdict_cmd), shell=True, capture_output=False, check=True)
        else:
            print ('\033[31m %s.dict  exists. Not creating. \033[0m' % genome_fasta)

    def haplotype_caller(self, genome_fasta, infile, outfile):
        cmd = [self.cmd, "HaplotypeCaller",
               "-R", genome_fasta, "-I", infile,
               "-ploidy", "1", "-stand-call-conf", "30",
               "-O", outfile]
        print("++++++ GATK HaplotypeCaller Comnand: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def select_variants(self, genome_fasta, infile, outfile, select_type):
        cmd = [self.cmd, "SelectVariants",
               "-R", genome_fasta, "-V", infile,
               "-select-type", select_type,
               "-O", outfile]
        print("++++++ Select Variants: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def variant_filtration(self, genome_fasta, infile, outfile,
                           filter_name, filter_expression):
        cmd = [self.cmd, "VariantFiltration",
               "-R", genome_fasta, "-V", infile,
               "--filter-expression", "\"%s\"" % filter_expression,
               "--filter-name", "\"%s\"" % filter_name,
               "-O", outfile]
        print("++++++ Applying filters: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
