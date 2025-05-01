"""varscan.py - Interface to the varscan command line tool

This module provides a Python interface to the VarScan command
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess


# Get file names universally used
def get_pileup_file(varscan_results, exp_name):
    return "%s.pileup" % get_varscan_files_path(varscan_results, exp_name)


def get_snps_file(varscan_results, exp_name):
    return '%s_varscan_snps_final.vcf' % get_varscan_files_path(varscan_results, exp_name)


def get_indels_file(varscan_results, exp_name):
    return '%s_varscan_inds_final.vcf' % get_varscan_files_path(varscan_results, exp_name)


def get_cns_file(varscan_results, exp_name):
    return '%s_varscan_cns_final.vcf' % get_varscan_files_path(varscan_results, exp_name)


def get_varscan_files_path(varscan_results, exp_name):
    """Get the path to the VarScan files"""
    return os.path.join(varscan_results, exp_name)


class VarScan:
    """This interface was created based on VarScan 2.4.0"""

    def __init__(self, config):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.config = config
        self.cmd = config['tools']['varscan']['cmd']

    def mpileup2snp(self, varscan_results, exp_name):
        """VarScan for SNPs"""
        pileup_file = get_pileup_file(varscan_results, exp_name)
        outfile = get_snps_file(varscan_results, exp_name)
        conf = self.config["tools"]["varscan"]["mpileup2snp"]

        cmd = [self.cmd, "mpileup2snp", pileup_file,
               "--output-vcf", "1",
               "--min-coverage", str(conf["min-coverage"]),
               "--min-reads2", str(conf["min-reads2"]),
               "--min-avg-qual", str(conf["min-avg-qual"]),
               "-min-freq-for-hom", str(conf["min-freq-for-hom"]),
               "--strand-filter", "1" if conf["strand-filter"] else "0",
               ">",
               outfile]
        command = ' '.join(cmd)
        print("++++++ Varscan for SNPs: '%s'" % command)
        compl_proc = subprocess.run(command, shell=True, capture_output=False,
                                    check=True)

    def mpileup2indel(self, varscan_results, exp_name):
        """VarScan for indels"""
        pileup_file = get_pileup_file(varscan_results, exp_name)
        outfile = get_indels_file(varscan_results, exp_name)
        conf = self.config["tools"]["varscan"]["mpileup2snp"]

        cmd = [self.cmd, "mpileup2indel", pileup_file,
               "--output-vcf", "1",
               "--min-coverage", str(conf["min-coverage"]),
               "--min-reads2", str(conf["min-reads2"]),
               "--min-avg-qual", str(conf["min-avg-qual"]),
               "-min-freq-for-hom", str(conf["min-freq-for-hom"]),
               "--strand-filter", "1" if conf["strand-filter"] else "0",
               ">",
               outfile]
        command = ' '.join(cmd)
        print("++++++ Varscan for INDELS: '%s'" % command)
        compl_proc = subprocess.run(command, shell=True, capture_output=False,
                                    check=True)

    def mpileup2cns(self, varscan_results, exp_name):
        """VarScan for everything"""
        pileup_file = get_pileup_file(varscan_results, exp_name)
        outfile = get_cns_file(varscan_results, exp_name)
        conf = self.config["tools"]["varscan"]["mpileup2snp"]

        cmd = [self.cmd, "mpileup2cns", pileup_file,
               "--output-vcf", "1",
               "--min-coverage", str(conf["min-coverage"]),
               "--min-reads2", str(conf["min-reads2"]),
               "--min-avg-qual", str(conf["min-avg-qual"]),
               "-min-freq-for-hom", str(conf["min-freq-for-hom"]),
               "--strand-filter", "1" if conf["strand-filter"] else "0",
               ">",
               outfile]
        command = ' '.join(cmd)
        print("++++++ Varscan for CNSs: '%s'" % command)
        compl_proc = subprocess.run(command, shell=True, capture_output=False,
                                    check=True)
