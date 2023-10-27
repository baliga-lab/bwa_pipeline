#!/bin/bash
# This script is the high level script to
# submit an entire directory of SRA samples to the cluster for
# BWA processing

# The steps are
# 0. for each sample in the input directory
#   1. Create FASTQ files using fastq-dump --split-3 <srafile> into a dedicated directory
#   2. run BWA pipleline on the data directory
#      Example
#      ./bwa_pipeline.py /proj/omics4tb2/wwu/BWA_pipeline/code/evan_data SRR11786873 newbwa_results_deleteme
#   3. delete the FASTQ files
SAMPLES=$(cat TODO_SAMPLEBATCHES/bwainput_samples-aa)
DATA_DIR="/proj/omics4tb2/wwu/BWA_HUGE_DATA"
FAILED_SAMPLES="/proj/omics4tb2/wwu/BWA_pipeline/evan_failed_samples.txt"
RESULT_DIR="/proj/omics4tb2/wwu/BWA_HUGE_RESULTS"

if [ ! -d "$DATA_DIR" ]; then
    echo "Creating $DATA_DIR"
    mkdir $DATA_DIR
fi
if [ ! -d "$RESULT_DIR" ]; then
    echo "Creating $RESULT_DIR"
    mkdir $RESULT_DIR
fi

for line in $SAMPLES
do
    # sample SRALite file
    SRA_FILE="/proj/omics4tb2/wwu/SRA_DOWNLOADS/SRALite_REAL/$line"
    echo $SRA_FILE
    arr=(${line//./ })
    samplename=${arr[0]}
    sampledir="$DATA_DIR/$samplename"
    mkdir -p $sampledir
    pushd $sampledir
    if fastq-dump --split-3 $SRA_FILE ; then
	batch_job="$samplename.sh"
	echo "./make_bwa_job.py $DATA_DIR $samplename $RESULT_DIR > $batch_job"
	popd
	./make_bwa_job.py $DATA_DIR $samplename $RESULT_DIR > $batch_job && sbatch $batch_job
    else
	echo "FAILURE !!!"
	echo $samplename >> $FAILED_SAMPLES
	popd
    fi
done
