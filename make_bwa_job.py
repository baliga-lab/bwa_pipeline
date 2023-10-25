#!/usr/bin/env python3
import jinja2
import argparse

TEMPLATE = """#!/bin/bash

#SBATCH -J bwa_pipeline
#SBATCH -o {{log_dir}}/"%j".out
#SBATCH -e {{log_dir}}/"%j".out

./bwa_pipeline.py --organism mtb {{datadir}} {{sample}} {{resultdir}} && rm -rf {{datadir}}/{{sample}}
"""

DESCRIPTION="""make_bwa_job - create a slurm job for BWA"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('datadir', help="data directory full path")
    parser.add_argument('sample', help="sample directory within data directory")
    parser.add_argument('resultdir', help="sample directory within data directory")
    args = parser.parse_args()
    templ = jinja2.Template(TEMPLATE)
    config = {}
    config['log_dir'] = "/proj/omics4tb2/wwu/BWA_pipeline/evan_slurm_log"
    config['datadir'] = args.datadir
    config['sample'] = args.sample
    config['resultdir'] = args.resultdir
    print(templ.render(config))
