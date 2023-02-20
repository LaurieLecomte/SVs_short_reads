#!/bin/bash

# Filter for SV calls tagged as PRECISE in smoove's calls

# valeria
# srun -c 1 -p ibis_small -J 03.5_smoove_filter -o log/03.5_smoove_filter_%j.log /bin/sh 01_scripts/03.5_smoove_filter.sh &

# manitou
# srun -c 1 -p small -J 03.5_smoove_filter -o log/03.5_smoove_filter_%j.log /bin/sh 01_scripts/03.5_smoove_filter.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"


# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Remove imprecise calls and BNDs
bcftools filter -i 'IMPRECISE=0 & SVTYPE!="BND"' $MERGED_DIR/smoove/smoove_merged_sorted.vcf.gz | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END | bcftools sort > $FILT_DIR/smoove/smoove_PRECISE.vcf
#tabix -p vcf $FILT_DIR/smoove/smoove_PRECISE.vcf.gz

#bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END $FILT_DIR/smoove/smoove_PRECISE.vcf > $FILT_DIR/smoove/smoove_PRECISE_simpl.vcf