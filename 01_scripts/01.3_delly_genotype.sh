#!/bin/bash

# 3rd step of SV calling with delly : genotype SV calls at each sites
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes : we have high coverage (16X)
# Following instructions for germline SV calling (https://github.com/dellytools/delly#germline-sv-calling)

# Activate the conda env before running: conda activate SVs_SR

# Manitou
# parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -c 1 -p small --time=1-00:00:00 --mem=20G -J 01.3_delly_genotype_{} -o log/01.3_delly_genotype_{}_%j.log /bin/sh 01_scripts/01.3_delly_genotype.sh {} &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

SAMPLE=$1
BAM="$BAM_DIR/"$SAMPLE".bam"

REGIONS_EX="02_infos/excl_chrs.txt"

# LOAD REQUIRED MODULES
#module load delly/1.1.6
#module load bcftools/1.13

# Create directory for genotyped calls
if [[ ! -d $CALLS_DIR/delly/geno ]]
then
  mkdir $CALLS_DIR/delly/geno
fi

# 1. Genotype SV sites across all samples
delly call -g $GENOME -v $CALLS_DIR/delly/merged_sites/merged_sites.bcf -o $CALLS_DIR/delly/geno/"$SAMPLE".geno.bcf $BAM
