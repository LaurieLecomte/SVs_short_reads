#!/bin/bash

# First step of SV calling with delly : call SV for each sample 
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes : we have high coverage (16X)
# Following instructions for germline SV calling (https://github.com/dellytools/delly#germline-sv-calling)

# parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -c 1 --mem=20G -p ibis_medium --time=7-00:00 -J 01.1_delly_call_{} -o log/01.1_delly_call_{}_%j.log /bin/sh 01_scripts/01.1_delly_call.sh {} &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

SAMPLE=$1

REGIONS_EX="02_infos/excl_chrs.txt"

# LOAD REQUIRED MODULES
module load delly/1.1.6
module load bcftools/1.15

# Create directory for raw calls
if [[ ! -d $CALLS_DIR/raw ]]
then
  mkdir $CALLS_DIR/raw
fi

# 1. Run delly on each file
delly call -g $GENOME $BAM_DIR/"$SAMPLE".bam -o $CALLS_DIR/delly/raw/"$SAMPLE".bcf -x $REGIONS_EX



