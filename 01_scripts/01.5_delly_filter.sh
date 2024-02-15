#!/bin/bash

# Filter for SV calls tagged as PASS and PRECISE in delly's calls

# valeria
# srun -c 1 -p ibis_small -J 01.5_delly_filter --mem=20G --time=1-00:00:00 -o log/01.5_delly_filter_%j.log /bin/sh 01_scripts/01.5_delly_filter.sh &

# manitou
# srun -c 1 -p small -J 01.5_delly_filter -o log/01.5_delly_filter_%j.log /bin/sh 01_scripts/01.5_delly_filter.sh &

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
#module load delly/1.1.6
#module load bcftools/1.13

# 1. Filter for PASS and PRECISE calls, remove BNDs, then extract required fields
bcftools filter -i 'FILTER="PASS" & PRECISE=1 & SVTYPE!="BND"' $MERGED_DIR/delly/delly_merged_sorted.vcf.gz | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END,INFO/CONSENSUS | bcftools sort > $FILT_DIR/delly/delly_PASS_PRECISE.vcf
#tabix -p vcf $FILT_DIR/delly/delly_PASS_PRECISE.vcf.gz

#bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END,INFO/CONSENSUS $FILT_DIR/delly/delly_PASS_PRECISE.vcf > $FILT_DIR/delly/delly_PASS_PRECISE_simpl.vcf