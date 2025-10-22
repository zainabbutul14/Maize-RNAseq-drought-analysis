#!/bin/bash
# ----------------------------------------------------------
# RNA-seq Analysis Pipeline for Maize Drought Project
# Author: Zainab Butul
# Description: Full bash pipeline from raw data download to quantification 
#All steps were executed manually through the terminal.
#Below are the exact commands used to perform the RNA-seq analysis.
# ----------------------------------------------------------

# Create project directory
mkdir maize_rnaseq_project
cd maize_rnaseq_project

# Make a folder to save FASTA files
mkdir -p ~/Users/zainabbutul/maize_rnaseq_project/ALL_FASTQ
cd ~/Users/zainabbutul/maize_rnaseq_project/ALL_FASTQ

#  Install essential tools 
brew install sratoolkit
which fasterq-dump

#  Download RNA-seq data
awk -F',' 'NR>1 {print $1}' Project-1_metadata.csv > run.txt
while read srr; do
    echo "Downloading $srr ..."
    fasterq-dump $srr -O ~/maize_rnaseq_project/ALL_FASTQ
    echo "$srr done"
done < run_ids.txt
ls ~/Maize_RNAseq_Project/ALL_FASTQ

# Quality Control & Trimming 
python3-m venv.env
source .venv/bin/activate
pip install --upgrade pip
pip install pyfastx cutadapt pandas
brew install fastqc fastp 
python3 qc_trimming.py

# Quantification
brew install kallisto
kallisto version
kallisto index -i /Users/zainabbutul/maize_rnaseq_project/reference/maize_index.idx/ /Users/zainabbutul/maize_rnaseq_project/reference/transcripts.fa

for fq in /Users/zainabbutul/maize_rnaseq_project/ALL_FASTQ/*.fastq
do
    sample=$(basename $fq .fastq)
    outdir=/Users/zainabbutul/maize_rnaseq_project/QUANT/$sample
    mkdir -p $outdir

    echo "Quantifying $sample ..."
    kallisto quant -i /Users/zainabbutul/maize_rnaseq_project/reference/maize_index.idx \
        -o $outdir --single -l 200 -s 20 $fq
done

# post-processing
python3 quantification.py
python3 gtf.py