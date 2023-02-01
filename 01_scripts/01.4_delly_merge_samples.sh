#!/bin/bash

# 4th step of SV calling with delly : merge sample vcf 
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes : we have high coverage (16X)
# Following instructions for germline SV calling (https://github.com/dellytools/delly#germline-sv-calling)

# srun -c 1 -p ibis_small -J 01.4_delly_merge_samples -o log/01.4_delly_merge_samples_%j.log /bin/sh 01_scripts/01.4_delly_merge_samples.sh &

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

# Remove previous list from previous trials
if [[ -f $BCF_GENO_LIST ]]
then
  rm $BCF_GENO_LIST
fi


# 1. Make a list of bcf files to merge
ls -1 $CALLS_DIR/delly/geno/*.bcf > $BCF_GENO_LIST

# 2. Merge all genotyped samples to get a single VCF/BCF using bcftools merge
bcftools merge -m id -O b -o $CALLS_DIR/delly/merged_samples.bcf -l $BCF_GENO_LIST

# 3. Convert bcf to vcf using bcftools, add tags and sort
bcftools view $CALLS_DIR/delly/merged_samples.bcf | bcftools +fill-tags | bcftools sort -Oz > $MERGED_DIR/delly/delly_merged_sorted.vcf.gz
