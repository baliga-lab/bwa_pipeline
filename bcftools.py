"""bcftools.py - Interface to the bcftools command line tools

This module provides a Python interface to various bcftools commands.
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class BcfTools:
    """This interface was created based on bcftools 1.18"""

    def __init__(self, bcftools_path, tabix_path, bgzip_path):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.bcftools = bcftools_path
        self.tabix = tabix_path
        self.bgzip = bgzip_path


    def variant_calling_mpileup(self, genome_fasta, infile, outfile, num_threads=16):
        cmd = [self.bcftools, "mpileup", "--threads", str(num_threads), "-Ou",
               "-f", genome_fasta, infile, "|",  # note the "|" operator !!!
               self.bcftools, "call", "-vmO", "v",
               "--ploidy", "1", "-o", outfile]
        print("++++++ Variant Calling mpileup: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)


    def filter_variants(self, infile, outfile,
                        output_type="v", soft_filter="LOWQUAL", include_expression="\"QUAL>10\""):
        cmd = [self.bcftools, "filter", "-O", output_type,
               "-o", outfile, "-s", soft_filter,
               "-i", include_expression, infile]
        print("++++++ Variant Calling filtering: ", ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def view(self, infile, outfile, output_type="v", allele_frequency=0.01):
        cmd = [self.bcftools, "view", "-O", output_type,
               "-o", outfile, "-q", str(allele_frequency),
               infile]
        print("++++++ FINAL view filtering: ", ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
