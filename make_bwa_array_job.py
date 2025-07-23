#!/usr/bin/env python3
import jinja2
import argparse
import json

TEMPLATE = """#!/bin/bash

#SBATCH -J bwa_pipeline2
#SBATCH -o {{log_dir}}/"%A"."%a".out
#SBATCH -e {{log_dir}}/"%A"."%a".out
#SBATCH --array={{array_range}}

data_folders=({{data_folders}})
data_folder=${data_folders[$SLURM_ARRAY_TASK_ID]}

./bwa_pipeline.py {{config}} $data_folder
"""

DESCRIPTION="""make_bwa_job - create a slurm job for BWA"""

ARRAY_MAX_TASKS = 10
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('config', help="configuration file")
    parser.add_argument('samplefile', help="file containing sample names")
    args = parser.parse_args()
    templ = jinja2.Template(TEMPLATE)
    with open(args.config) as infile:
        config = json.load(infile)
    config["config"] = args.config
    samples = []
    with open(args.samplefile) as infile:
        for line in infile:
            samples.append(line.strip())
    config['data_folders'] = " ".join(["\"%s\"" % x for x in samples])

    array_max_task_spec = "%%%d" % ARRAY_MAX_TASKS
    config["array_range"] = "0-%d%s" % (len(samples) - 1, array_max_task_spec)
    print(templ.render(config))
