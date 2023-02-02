#!/bin/bash

# 4th step of SV calling with smoove : merge sites across samples 
# srun -c 1 --mem=20G -p small -J 04_paste -o log/04_paste_%j.log /bin/sh 01_scripts/04_paste.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

# LOAD REQUIRED MODULES
module load bcftools

# 1. Paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.
smoove paste --outdir $MERGED_DIR/smoove --name merged $CALLS_DIR/smoove/geno/*.vcf.gz

# 2. Add tags and sort
bcftools view $MERGED_DIR/smoove/merged.smoove.square.vcf.gz | bcftools +fill-tags | bcftools sort -Oz > $MERGED_DIR/smoove/smoove_merged_sorted.vcf.gz

