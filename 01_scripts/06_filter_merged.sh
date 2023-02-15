#!/bin/bash

# Filter merged output to keep calls supported by at least 2 callers, add tags for further filtering if needed

# srun -c 2 -p small -J 06_filter_merged -o log/06_filter_merged_%j.log /bin/sh 01_scripts/06_filter_merged.sh &
 
# VARIABLES
GENOME="03_genome/genome.fasta"
CALLS_DIR="05_calls"
FILT_DIR="07_filtered"

MERGED_UNION_DIR="08_merged_union"
FILT_UNION_DIR="09_filtered_union"

DELLY_VCF="$FILT_DIR/delly/delly_PASS_PRECISE.vcf"
MANTA_VCF="$FILT_DIR/manta/manta_PASS_PRECISE.vcf"
SMOOVE_VCF="$FILT_DIR/smoove/smoove_PRECISE.vcf"

VCF_LIST="02_infos/callers_VCFs.txt"
MERGED_VCF="$MERGED_UNION_DIR/merged_delly_manta_smoove.vcf"

REGIONS_EX="02_infos/excl_chrs.txt"

CPU=2

# LOAD REQUIRED MODULES
module load bcftools

# 1. Filter for SVs over 50 bp and supported by at least 2 tools
bcftools filter -i 'SUPP!="1" & ABS(SVLEN) >= 50' $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".sorted.vcf --threads $CPU | bcftools +fill-tags --threads $CPU | bcftools sort > $FILT_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)"_SUPP2.vcf

# 2. Extract indels (<50 bp) supported by at least 2 tools, for reference
bcftools filter -i 'SUPP!="1" & ABS(SVLEN) < 50' $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".sorted.vcf --threads $CPU | bcftools +fill-tags --threads $CPU | bcftools sort > $FILT_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)"_SUPP2_indels.vcf