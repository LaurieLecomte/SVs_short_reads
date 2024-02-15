#!/bin/bash

# Filter for SV calls tagged as PASS and PRECISE in manta's calls

# valeria
# srun -c 1 -p ibis_small --mem=20G --time=1-00:00:00 -J 02.3_manta_filter -o log/02.3_manta_filter_%j.log /bin/sh 01_scripts/02.3_manta_filter.sh &

# manitou
# srun -c 1 -p small -J 02.3_manta_filter -o log/02.3_manta_filter_%j.log /bin/sh 01_scripts/02.3_manta_filter.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"


# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Filter for PASS and PRECISE calls, remove BNDs if any, then extract required fields
bcftools filter -i 'FILTER="PASS" & IMPRECISE=0 & SVTYPE!="BND"' $MERGED_DIR/manta/manta_merged_sorted.vcf.gz | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END | bcftools sort > $FILT_DIR/manta/manta_PASS_PRECISE.vcf
#tabix -p vcf $FILT_DIR/manta/manta_PASS_PRECISE.vcf.gz

#bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END $FILT_DIR/manta/manta_PASS_PRECISE.vcf > $FILT_DIR/manta/manta_PASS_PRECISE_simpl.vcf