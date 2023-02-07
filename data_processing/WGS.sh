#!/bin/bash

BAM_DIRECTORY=$1
OUT_DIRECTORY=$2
MITOCALC_DIRECTORY=$3

## 1. Telomere length (assuming conda is already installed)
### install telseq
conda install -c bioconda telseq
### create bam file list
ls $BAM_DIRECTORY/*.bam > $BAM_DIRECTORY/telseq_bamlist.txt
### run telseq (one file for all samples)
telseq -f $BAM_DIRECTORY/telseq_bamlist.txt -o $OUT_DIRECTORY/WGS/telseq_results.txt

## 2. Mitochondrial DNA copy number
cd $MITOCALC_DIRECTORY
### download and unzip the files
wget https://github.com/HSGU-NIA/mitoAnalyzer/raw/master/fastMitoCalc.zip
unzip fastMitoCalc.zip
### run fastmitocalc (one file for each sample)
for s in $BAM_DIRECTORY/*.bam 
do
perl $MITOCALC_DIRECTORY/fastMitoCalc.pl -e chr -m chrM -f $BAM_DIRECTORY/"$s".bam -w $OUT_DIRECTORY -p $MITOCALC_DIRECTORY/BaseCoverage
done
