#!/bin/bash

# This script produces required bed files beforehand to indicate regions where SVs must NOT or MUST be called
# If some chromosomes need to be removed from SV calling, they should be added to 02_infos/excl_chrs.txt
# WARNING : the 02_infos/excl_chrs.txt file mush be encoded in linux format, otherwise grep won't grep, AND have a newline at the end

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"

REGIONS_EX="02_infos/excl_chrs.txt"

# LOAD REQUIRED MODULES
module load htslib
module load samtools

if [[ -f 02_infos/excl_chrs.bed ]]
then
  rm 02_infos/excl_chrs.bed
fi
if [[ -f 02_infos/chrs.bed ]]
then
  rm 02_infos/chrs.bed
fi


# 1. Create bed from genome index
less "$GENOME".fai | awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' > "$GENOME".bed

# 2. Generate bed for excluded chromosomes or contigs (for smoove)
less $REGIONS_EX | while read REGION || [ -n "$line" ]
do
  grep -Fw $REGION "$GENOME".bed >> 02_infos/excl_chrs.bed
done

bgzip -f 02_infos/excl_chrs.bed
tabix -f -p bed 02_infos/excl_chrs.bed.gz


# 3. Generate bed for chromosomes or contigs for which SVs are to be called  (for manta)
less $CHR_LIST | while read CHR || [ -n "$line" ]
do
  grep -Fw $CHR "$GENOME".bed >> 02_infos/chrs.bed
done

bgzip -f 02_infos/chrs.bed
tabix -f -p bed 02_infos/chrs.bed.gz