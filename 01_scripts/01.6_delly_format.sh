#!/bin/bash

# Format SV calls by delly prior to merging across callers

# srun -c 1 -p ibis_small -J 01.6_delly_format -o log/01.6_delly_format_%j.log /bin/sh 01_scripts/01.6_delly_format.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"




# 1. Add explicit ALT sequence for INVs, DELs


# 2. Add SVLEN field


# 3. Simplify VCF
