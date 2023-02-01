#!/bin/bash

# Call SV in all samples, parallelized by chr

# parallel -a 02_infos/chr_list.txt -k -j 10 srun -c 4 --mem=20G -p medium --time=7-00:00 -J 01_call_{} -log/01_call_{}_%j.log /bin/sh 01_scripts/01_call.sh {} &
# srun -c 8 -p ibis_medium --time=7-00:00 -J 02.1_manta_call -o log/02.1_manta_call_%j.log /bin/sh 01_scripts/02.1_manta_call.sh  &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

REGIONS="02_infos/chrs.bed.gz"

BAM_LIST=$(for file in $(ls $BAM_DIR/*.bam); do echo '--bam' "$file" ; done)

CPU=8

# LOAD REQUIRED MODULES
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
module load manta/1.6.0
module load bcftools

## Paths and exec locations for running manta
SAMTOOLS_PATH=$(which samtools)

# Increase opened file number limit
ulimit -S -n 2048


# 1. Workflow configuration : set reference genome and samples for which SV are to be called in order to generate an executable (runWorkflow.py)
configManta.py --referenceFasta $GENOME --runDir $CALLS_DIR/manta --callRegions $REGIONS $(echo $BAM_LIST)

# 3. Launch resulting executable
## -j controls the number of cores/nodes
$CALLS_DIR/manta/runWorkflow.py -j $CPU

# 4. Convert BNDs to INVs : the convertInversion.py script changes BNDs to INVs and adds a SVLEN field to these SVs
01_scripts/utils/convertInversion.py $SAMTOOLS_PATH $GENOME $CALLS_DIR/manta/results/variants/diploidSV.vcf.gz | bgzip > $CALLS_DIR/manta/results/variants/diploidSV_converted.vcf.gz

# 5. Sort and rename output
bcftools sort $CALLS_DIR/manta/results/variants/diploidSV_converted.vcf.gz -Oz > $CALLS_DIR/manta/manta_sorted.vcf.gz