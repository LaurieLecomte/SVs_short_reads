#!/bin/bash

# Merge SV calls across the 3 callers used
# Jasmine must be installed in current env or session

# As of 20230214, can be done on manitou only, because Jasmine has not been installed yet on valeria
# valeria
# srun -c 2 -p ibis_small -J 04_merge -o log/04_merge_%j.log /bin/sh 01_scripts/04_merge_callers.sh &

# manitou
# srun -c 2 -p small -J 04_merge -o log/04_merge_%j.log /bin/sh 01_scripts/04_merge_callers.sh &
 
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

CPU=2

# LOAD REQUIRED MODULES


# 0. Remove VCF list from previous trials if any
if [[ -f $VCF_LIST ]]
then
  rm $VCF_LIST
fi

# 1. List all merged VCFs
touch $VCF_LIST
 

echo $DELLY_VCF > $VCF_LIST
echo $MANTA_VCF >> $VCF_LIST
echo $SMOOVE_VCF >> $VCF_LIST

# 2. Merge SV calls accross samples, using predefined parameters 
jasmine file_list=$VCF_LIST out_file=$MERGED_VCF out_dir=$MERGED_UNION_DIR genome_file=$GENOME --ignore_strand --ignore_merged_inputs --normalize_type --output_genotypes --allow_intrasample --mutual_distance --max_dist_linear=0.25 --threads=$CPU
 

