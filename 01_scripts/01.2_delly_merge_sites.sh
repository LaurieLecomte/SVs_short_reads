#!/bin/bash

# 2nd step of SV calling with delly : merge sites 
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes : we have high coverage (16X)
# Following instructions for germline SV calling (https://github.com/dellytools/delly#germline-sv-calling)

# srun -p ibis_small --time=1-00:00:00 -J 01.2_delly_merge_sites -o log/01.2_delly_merge_sites_%j.log /bin/sh 01_scripts/01.2_delly_merge_sites.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS_EX="02_infos/excl_chrs.txt"

BCF_LIST="02_infos/bcf_list.txt"

# LOAD REQUIRED MODULES
#module load delly/1.1.6
#module load bcftools


# 0. Create directory for merged sites
if [[ ! -d $CALLS_DIR/delly/merged_sites ]]
then
  mkdir $CALLS_DIR/delly/merged_sites
fi

# 1. Make a list of bcf files to merge
## Remove previous list from previous trials
if [[ -f $BCF_LIST ]]
then
  rm $BCF_LIST
fi

ls -1 $CALLS_DIR/delly/raw/*.bcf > $BCF_LIST

# 2. Merge SV sites into a unified site list
delly merge -o $CALLS_DIR/delly/merged_sites/merged_sites.bcf $BCF_LIST