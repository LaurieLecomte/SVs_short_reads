#!/bin/bash

# 2nd step of SV calling with smoove : merge sites across samples

# RUN ON MANITOU ONLY
# srun -c 1 --mem=20G -p small -J 03.2_smoove_merge -o log/03.2_smoove_merge_%j.log /bin/sh 01_scripts/03.2_smoove_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

# LOAD REQUIRED MODULES
module load python/2.7 gsort/0.1.4 lumpy-sv/0.3.1 svtyper/0.7.1 smoove/0.2.7 
module load samtools/1.12
module load bcftools/1.13

# 0. Create directory for merged sites
if [[ ! -d $CALLS_DIR/smoove/merged_sites ]]
then
  mkdir $CALLS_DIR/smoove/merged_sites
fi

# 1. Get the union of sites across all samples 
smoove merge --name merged -f $GENOME --outdir $CALLS_DIR/smoove/merged_sites $CALLS_DIR/smoove/raw/*raw.vcf.gz