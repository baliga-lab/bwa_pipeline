"""varscan.py - Interface to the varscan command line tool

This module provides a Python interface to the VarScan command
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class VarScan:
    """This interface was created based on VarScan 2.4.0"""

    def __init__(self, cmd_path):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.cmd = cmd_path

    def mpileup2snp(self, infile, outfile):
        """VarScan for SNPs"""
        #cmd2 = '%s mpileup2snp %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_snps_final.vcf' % (config['tools']['varscan'],
        #                                                                                                                                                varscan_files_path,
        #                                                                                                                                                varscan_files_path)
        cmd = [self.cmd, "mpileup2snp", infile,
               "--output-vcf", "1", "--min-coverage", "8",
               "--min-reads2", "2", "--min-avg-qual", "30",
               "--strand-filter", "0", ">",
               outfile]
        print("++++++ Varscan for SNPs: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def mpileup2indel(self, infile, outfile):
        """VarScan for indels"""
        # cmd3 = '%s mpileup2indel %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_inds_final.vcf' % (config['tools']['varscan'],
        #                                                                                                                                                  varscan_files_path,
        #                                                                                                                                                  varscan_files_path)
        cmd = [self.cmd, "mpileup2indel", infile,
               "--output-vcf", "1", "--min-coverage", "8",
               "--min-reads2", "2", "--min-avg-qual", "30",
               "--strand-filter", "0", ">",
               outfile]
        print("++++++ Varscan for INDELS: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
