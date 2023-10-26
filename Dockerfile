# syntax=docker/dockerfile:1
FROM continuumio/anaconda3:latest

RUN apt update
RUN apt-get install -y git wget libncurses-dev libbz2-dev liblzma-dev unzip vim bzip2
RUN apt-get install -y gcc g++ zlib1g-dev build-essential
RUN apt-get install -y openjdk-11-jdk
RUN apt-get install -y r-base

WORKDIR /

# build SAMTools
RUN wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
RUN tar xvjf samtools-1.6.tar.bz2
WORKDIR /samtools-1.6
RUN ./configure
RUN make && make install
WORKDIR /

# build GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
RUN unzip gatk-4.3.0.0.zip
WORKDIR /gatk-4.3.0.0
RUN sed -e 's/python/python3/g' < gatk > gatk4
RUN mv gatk4 gatk
RUN chmod a+x gatk

# build BCFTools, this includes htslib, tabix and bgzip
WORKDIR /
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
RUN tar xvjf bcftools-1.18.tar.bz2
WORKDIR /bcftools-1.18
RUN ./configure
RUN make && make install
WORKDIR /bcftools-1.18/htslib-1.18
RUN make && make install

# BWA 0.7.17
WORKDIR /
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /bwa
RUN make

# snpEff-5.1d, unpacks to "snpEff"
WORKDIR /
RUN wget https://networks.systemsbiology.net/downloads/snpEff-5.1d.tar.gz
RUN tar xfz snpEff-5.1d.tar.gz

# Install jar files
RUN mkdir -p /jarfiles
WORKDIR /jarfiles
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O picard-2.27.5.jar
RUN wget https://github.com/dkoboldt/varscan/releases/download/2.4.0/VarScan.v2.4.0.jar

# Python dependencies
RUN pip install globalsearch

# Install trimgalore and fastqc
RUN pip install cutadapt
RUN conda install -c bioconda fastqc

# Boto3 for AWS integration
RUN conda install -y boto3

WORKDIR /
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -O trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz

# adjust path
RUN echo "PATH=\$PATH:/gatk-4.3.0.0:/bwa:/TrimGalore-0.6.10" >> ~/.bashrc
WORKDIR /

RUN mkdir -p /bwa_data
RUN mkdir -p /bwa_results
RUN mkdir -p /bwa_tmp

# Install pipeline software and data
RUN git clone https://github.com/baliga-lab/bwa_pipeline.git
WORKDIR /bwa_pipeline
RUN wget https://networks.systemsbiology.net/downloads/bwa_mtb_reference_genome-20231026.tar.gz
RUN tar xvf bwa_mtb_reference_genome-20231026.tar.gz
RUN rm bwa_mtb_reference_genome-20231026.tar.gz
