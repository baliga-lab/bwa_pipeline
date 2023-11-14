#!/bin/bash
# This script is the high level script to download and run data
# for a given sample

# The steps are
# 0. download the SRA file for the sample
# 1. create FASTQ files for the sample using fastq-dump --split-3
# 2. run BWA pipleline on the data directory

sra_url=$(./get_sra_url.py mtb_downloads_sralite.csv $1 2>&1)
if [ $? -eq "0" ] ; then
    target_dir=/bwa_data/$1
    mkdir -p $target_dir
    # clean up the target dir if necessary
    rm $target_dir/*
    echo "downloading $sra_url..."
    wget --directory-prefix=$target_dir/ $sra_url
    fasterq-dump -f -O $target_dir --split-3 $target_dir/*
    ./bwa_pipeline.py --organism mtb --config bwa_config.json $1
else
    echo "could not download $sra_url"
fi
