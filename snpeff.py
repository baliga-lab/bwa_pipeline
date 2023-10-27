"""snpEff.py - Interface to the snpEff command line tool

This module provides a Python interface to the snpEff command
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class SnpEff:
    """This interface was created based on snpEff 5.1d"""

    def __init__(self, cmd_path):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.cmd = cmd_path

    def format_eff(self, stats_file, snpeff_db, infile, outfile,
                   output_format="gatk"):
        #cmd1 ='%s -ud 0 -classic -csvStats %s -geneId -lof -v -formatEff -o gatk %s %s > %s' % (config['tools']['snpeff'],
        #                                                                                        snpeff_stats, snpeff_db, vcf_file_renamed, snpeff_vcf)
        cmd = [self.cmd, "-ud", "0", "-classic",
               "-csvStats", stats_file,
               "-geneId", "-lof", "-v",
               "-formatEff", "-o", "gatk",
               snpeff_db, infile, ">", outfile]
        print("++++++ Running SNPEff Formateff command: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
