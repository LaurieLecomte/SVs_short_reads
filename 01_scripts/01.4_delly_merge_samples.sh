#!/bin/bash

# 4th step of SV calling with delly : merge sample vcf 
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes : we have high coverage (16X)
# Following instructions for germline SV calling (https://github.com/dellytools/delly#germline-sv-calling)

# srun -c 1 -p ibis_small --mem=20G --time=1-00:00:00 -J 01.4_delly_merge_samples -o log/01.4_delly_merge_samples_%j.log 01_scripts/01.4_delly_merge_samples.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS_EX="02_infos/excl_chrs.txt"

BCF_GENO_LIST="02_infos/delly_bcf_list.txt"

# LOAD REQUIRED MODULES
#module load delly/1.1.6
#module load bcftools


# 1. Make a list of bcf files to merge
## Remove previous list from previous trials
if [[ -f $BCF_GENO_LIST ]]
then
  rm $BCF_GENO_LIST
fi

ls -1 $CALLS_DIR/delly/geno/*.bcf > $BCF_GENO_LIST

# 2. Merge all genotyped samples to get a single VCF/BCF using bcftools merge
bcftools merge -m id -O b -o $MERGED_DIR/delly/merged_samples.bcf -l $BCF_GENO_LIST

# 3. Convert bcf to vcf using bcftools, add tags, sort and remove excluded contigs from header
bcftools view $MERGED_DIR/delly/merged_samples.bcf | grep -vFf $REGIONS_EX  | bcftools sort -Oz > $MERGED_DIR/delly/delly_merged_sorted.vcf.gz

# Clean up
rm $MERGED_DIR/delly/merged_samples.bcf