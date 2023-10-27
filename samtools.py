"""samtools.py - Interface to the samtools command line tools

This module provides a Python interface to various samtools commands.
The intention is to facilitate parameterization and compatibility
issues between the various versions.
"""
import os
import subprocess

class SamTools:
    """This interface was created based on samtools 1.6"""

    def __init__(self, cmd_path):
        """Constructor.
        @param cmd_path path to the samtools command
        """
        self.cmd = cmd_path

    def faidx(self, genome_fasta):
        if not os.path.exists('%s.fai' % genome_fasta):
            faidx_cmd = [self.cmd, 'faidx', genome_fasta]
            print ("\033[31m '%s.fai' DOES NOT exist, creating. \033[0m" % genome_fasta)
            compl_proc = subprocess.run(' '.join(faidx_cmd), shell=True, capture_output=False, check=True)
        else:
            print ("\033[31m '%s.fai' exists. Not creating. \033[0m" % genome_fasta)

    def index(self, infile):
        """run the samtools index command"""
        index_cmd = [self.cmd, "index", infile]
        print ("++++++ Samtools Index Command: '%s'" % " ".join(index_cmd))
        compl_proc = subprocess.run(' '.join(index_cmd), shell=True, capture_output=False, check=True)

    def fixmate(self, infile, outfile, out_format="bam"):
        fixmate_cmd = [self.cmd, "fixmate",
                       "-O", out_format, infile, outfile]
        print ("++++++ Samtools Fixmate Command: '%s'" % ' '.join(fixmate_cmd))
        compl_proc = subprocess.run(' '.join(fixmate_cmd), shell=True, capture_output=False, check=True)

    def sort(self, infile, outfile, tmp_prefix, num_threads=8, out_format="bam"):
        sort_cmd = [self.cmd, "sort",
                    "-@", str(num_threads), "-O", out_format,
                    "-o", outfile, "-T", tmp_prefix, infile]
        print ("++++++ Samtools Sort Command: '%s'" % " ".join(sort_cmd))
        compl_proc = subprocess.run(' '.join(sort_cmd), shell=True, capture_output=False, check=True)

    def mpileup(self, genome_fasta, infile, outfile,
                no_baq=True):
        # original, but the --input-fmt-option nthreads=8 is weird
        #cmd1 = '%s mpileup --input-fmt-option nthreads=8 -B -f %s -o %s.pileup %s_marked.bam' % (config['tools']['samtools'],
        #                                                                                         genome_fasta, varscan_files_path, alignment_files_path)
        cmd = [self.cmd, "mpileup"]
        if no_baq:
            cmd.append("-B")
        cmd.extend(["-f", genome_fasta, "-o", outfile, infile])
        print("++++++ Samtools Mpileup: '%s'" % ' '.join(cmd))
        compl_proc = subprocess.run(' '.join(cmd), shell=True, capture_output=False, check=True)
