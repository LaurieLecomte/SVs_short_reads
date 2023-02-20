#!/bin/bash

# Call SV in all samples, parallelized by chr

# parallel -a 02_infos/chr_list.txt -k -j 10 srun -c 4 --mem=20G -p ibis_small --time=1-00:00 -J 02.1_manta_call_{} -o log/02.1_manta_call_{}_%j.log /bin/sh 01_scripts/02.1_manta_call.sh {} &
# srun -c 8 -p ibis_medium --time=7-00:00 -J 02.1_manta_call -o log/02.1_manta_call_%j.log /bin/sh 01_scripts/02.1_manta_call.sh  &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS="02_infos/chrs.bed.gz"

CHR=$1

BAM_LIST=$(for file in $(ls $BAM_DIR/*.bam); do echo '--bam' "$file" ; done)


CPU=4

# LOAD REQUIRED MODULES
#module load gcc python/3.10 manta/1.6.0
#module load bcftools/1.13
#module load samtools/1.15

## Paths and exec locations for running manta
SAMTOOLS_PATH=$(which samtools)

# Increase opened file number limit
ulimit -S -n 2048

# 0. Create output directory
if [[ ! -d $CALLS_DIR/manta/"$CHR" ]]
then
  mkdir $CALLS_DIR/manta/"$CHR"
else
  rm -r $CALLS_DIR/manta/"$CHR"/*
fi


# 1. Create bed for given chromosome from all chromosomes bed file
zless $REGIONS | grep -Fw "$CHR" > 02_infos/"$CHR".bed
bgzip 02_infos/"$CHR".bed -f
tabix -p bed 02_infos/"$CHR".bed.gz -f

# 1. Workflow configuration : set reference genome and samples for which SV are to be called in order to generate an executable (runWorkflow.py)
#configManta.py --referenceFasta $GENOME --runDir $CALLS_DIR/manta/$CHR --callRegions 02_infos/"$CHR".bed.gz $(echo $BAM_LIST)
configManta.py --referenceFasta $GENOME --runDir $CALLS_DIR/manta/$CHR --callRegions 02_infos/"$CHR".bed.gz  $BAM_LIST

# 3. Launch resulting executable
## -j controls the number of cores/nodes
$CALLS_DIR/manta/$CHR/runWorkflow.py -j $CPU

# 4. Convert BNDs to INVs : the convertInversion.py script changes BNDs to INVs and adds a SVLEN field to these SVs
01_scripts/utils/convertInversion.py $SAMTOOLS_PATH $GENOME $CALLS_DIR/manta/$CHR/results/variants/diploidSV.vcf.gz | bgzip > $CALLS_DIR/manta/$CHR/results/variants/diploidSV_converted.vcf.gz

# 5. Sort and rename output
bcftools sort $CALLS_DIR/manta/$CHR/results/variants/diploidSV_converted.vcf.gz -Oz > $CALLS_DIR/manta/manta_sorted_"$CHR".vcf.gz

# Clean up 
rm 02_infos/"$CHR".bed.gz
rm 02_infos/"$CHR".bed.gz.tbi