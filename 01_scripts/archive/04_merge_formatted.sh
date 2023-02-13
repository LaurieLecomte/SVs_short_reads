#!/bin/bash


 # srun -c 2 -p ibis_small -J 04_merge -o log/04_merge_%j.log /bin/sh 01_scripts/04_merge.sh &
 
# VARIABLES
GENOME="03_genome/genome.fasta"
CALLS_DIR="05_calls"
FILT_DIR="07_filtered"

MERGED_UNION_DIR="08_merged_union"
FILT_UNION_DIR="09_filtered_union"

DELLY_VCF="$FILT_DIR/delly/delly_PASS_PRECISE.vcf.gz"
MANTA_VCF="$FILT_DIR/manta/manta_PASS_PRECISE.vcf.gz"
SMOOVE_VCF="$FILT_DIR/smoove/smoove_PRECISE.vcf.gz"

VCF_LIST="02_infos/callers_VCFs.txt"

CPU=2
# LOAD REQUIRED MODULES



# 0. Remove VCF list from previous trials 
if [[ -f $VCF_LIST ]]
then
  rm $VCF_LIST
fi

#0. Initialize list of files to merge
touch $VCF_LIST

# 1. List all merged VCFs 
#echo "${DELLY_VCF%%.*}"_formatted.vcf > $VCF_LIST
#echo "${MANTA_VCF%%.*}"_formatted.vcf >> $VCF_LIST
#echo "${SMOOVE_VCF%%.*}"_formatted.vcf >> $VCF_LIST
echo "${MANTA_VCF%%.*}"_formatted.vcf > $VCF_LIST
echo "${DELLY_VCF%%.*}"_formatted.vcf >> $VCF_LIST
echo "${SMOOVE_VCF%%.*}"_formatted.vcf >> $VCF_LIST
 
# 2. Merge accross samples
#--max_dist_linear=0.1 --min_dist=50
jasmine file_list=$VCF_LIST out_file="$MERGED_UNION_DIR/manta_delly_smoove_distlin0.1_mutdist.vcf_mindist50" genome_file=$GENOME --ignore_strand --ignore_merged_inputs --normalize_type --output_genotypes --allow_intrasample --mutual_distance --max_dist_linear=0.1 --min_dist=50 --threads=$CPU