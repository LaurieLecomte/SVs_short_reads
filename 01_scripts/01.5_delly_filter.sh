#!/bin/bash

# Filter for SV calls tagged as PASS and PRECISE in delly's calls

# srun -c 1 -p ibis_small -J 01.5_delly_filter -o log/01.5_delly_filter_%j.log /bin/sh 01_scripts/01.5_delly_filter.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS_EX="02_infos/excl_chrs.txt"

BCF_GENO_LIST="02_infos/bcf_geno_list.txt"

# LOAD REQUIRED MODULES
module load delly/1.1.6
module load bcftools/1.15

# 1. Filter calls
bcftools filter -i 'FILTER="PASS" & PRECISE=1' $MERGED_DIR/delly/delly_merged_sorted.vcf.gz -Oz > $FILT_DIR/delly/delly_PASS_PRECISE.vcf.gz
tabix -p vcf $FILT_DIR/delly/delly_PASS_PRECISE.vcf.gz