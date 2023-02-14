#!/bin/bash

# Concat VCF from each chromosome
# srun -c 1 -p small -J 02.2_manta_merge -o log/02.2_manta_merge_%j.log /bin/sh 01_scripts/02.2_manta_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS="02_infos/chrs.bed.gz"

BAM_LIST=$(for file in $(ls $BAM_DIR/*.bam); do echo '--bam' "$file" ; done)

VCF_LIST="02_infos/manta_vcf_list.txt"

CPU=4

# LOAD REQUIRED MODULES
module load bcftools

# 1. Make a list of VCFs to concatenate
## Remove previous list from previous trials
if [[ -f $VCF_LIST ]]
then
  rm $VCF_LIST
fi

ls -1 $CALLS_DIR/manta/manta_sorted_* > $VCF_LIST

# 1. Concat, add tags and sort
bcftools concat -f $VCF_LIST  | bcftools sort -Oz > $MERGED_DIR/manta/manta_merged_sorted.vcf.gz