#!/bin/bash

# First step of SV calling with smoove : call SV for each sample 

# RUN ON MANITOU ONLY

# parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -p small -c 1 --mem=20G --time=1-00:00 -J 03.1_smoove_call_{} -o log/03.1_smoove_call_{}_%j.log /bin/sh 01_scripts/03.1_smoove_call.sh {} &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

SAMPLE=$1

# LOAD REQUIRED MODULES
module load python/2.7 gsort/0.1.4 samtools/1.12 lumpy-sv/0.3.1 svtyper/0.7.1 smoove/0.2.8


# Create directory for raw calls
if [[ ! -d $CALLS_DIR/smoove/raw ]]
then
  mkdir $CALLS_DIR/smoove/raw
fi

# 1. Call genotypes by parallelizing by SAMPLE
smoove call --outdir $CALLS_DIR/smoove/raw --name $SAMPLE --fasta $GENOME -p 1 --genotype $BAM_DIR/"$SAMPLE".bam --exclude 02_infos/excl_chrs.bed.gz

# 2. Rename for easier scripting
mv "$CALLS_DIR/smoove/raw/"$SAMPLE"-smoove.genotyped.vcf.gz" "$CALLS_DIR/smoove/raw/"$SAMPLE"_raw.vcf.gz"
mv "$CALLS_DIR/smoove/raw/"$SAMPLE"-smoove.genotyped.vcf.gz.csi" "$CALLS_DIR/smoove/raw/"$SAMPLE"_raw.vcf.gz.csi"