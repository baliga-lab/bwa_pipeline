"""snpEff.py - Interface to the snpEff command line tool

This module provides a Python interface to the snpEff command
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class SnpEff:
    """This interface was created based on snpEff 5.1d"""

    def __init__(self, snpeff, snpsift, vceff_opl):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.snpeff = snpeff
        self.snpsift = snpsift
        self.vceff_opl = vceff_opl

    def format_eff(self, stats_file, snpeff_db, infile, outfile,
                   output_format="gatk"):
        cmd = [self.snpeff, "-ud", "0", "-classic",
               "-csvStats", stats_file,
               "-geneId", "-lof", "-v",
               "-formatEff", "-o", "gatk",
               snpeff_db, infile, ">", outfile]
        print("++++++ Running SNPEff Formateff command: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def snpsift_filter(self, infile, outfile, filter_params):
        cmd = ["cat", infile, "|", self.snpsift, "filter",
               "-p", "\"%s\"" % filter_params,
               ">", outfile]
        print("++++++ Running SNPEff Filtering: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)

    def vceff_opl_extract_fields(self, infile, outfile, fields):
        """
        executes vcfeff_opl | snpeff extractFields
        The fields parameter is very flexible and needs to be properly double
        quoted to avoid shell expansion
        """
        cmd = ["cat", infile, "|", self.vceff_opl, "|", self.snpsift,
               "extractFields", "-",  # snpsift takes input from from stdin
               " ".join(fields),
               ">", outfile]
        print( "++++++ Running SNPEff Oneline final formatter: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
