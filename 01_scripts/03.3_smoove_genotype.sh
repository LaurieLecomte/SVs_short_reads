#!/bin/bash

# 3rd step of SV calling with smoove : genotype for all samples and all sites 

# parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -c 1 --mem=20G -p small -J 03.3_smoove_genotype_{} -o log/03.3_smoove_genotype_{}_%j.log /bin/sh 01_scripts/03.3_smoove_genotype.sh {} &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

SAMPLE=$1

# 0. Create directory for genotyped calls
if [[ ! -d $CALLS_DIR/smoove/geno ]]
then
  mkdir $CALLS_DIR/smoove/geno
fi


# 1. Genotype each sample at known SV sites from merged.sites.vcf.gz
smoove genotype -d -x -p 1 --name "$SAMPLE" --outdir $CALLS_DIR/smoove/geno --fasta $GENOME --vcf $CALLS_DIR/smoove/merged_sites/merged.sites.vcf.gz $BAM_DIR/"$SAMPLE".bam 