#!/bin/bash

# Filter for SV calls tagged as PASS and PRECISE in manta's calls

# srun -c 1 -p ibis_small -J 02.3_manta_filter -o log/02.3_manta_filter_%j.log /bin/sh 01_scripts/02.3_manta_filter.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"


# LOAD REQUIRED MODULES
module load bcftools/1.15

# 1. Filter calls
bcftools filter -i 'FILTER="PASS" & IMPRECISE=0' $MERGED_DIR/manta/manta_merged_sorted.vcf.gz | bcftools sort > $FILT_DIR/manta/manta_PASS_PRECISE.vcf
#tabix -p vcf $FILT_DIR/manta/manta_PASS_PRECISE.vcf.gz