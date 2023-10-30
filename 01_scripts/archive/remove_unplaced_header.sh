#!/bin/bash

# Removed unplaced scaffolds from VCF header

# VARIABLES
#INPUT_VCF=$1
#OUTPUT_VCF=$2

REGIONS_EX="02_infos/excl_chrs.txt"


##contig=<ID=

if [[ -f 02_infos/unplaced_contigs_lines.txt ]]
then
  rm 02_infos/unplaced_contigs_lines.txt
fi


# Make a file of strings to remove
less $REGIONS_EX | while read REGION
do
  echo "##contig=<ID="$REGION >> 02_infos/unplaced_contigs_lines.txt
done

# Remove these lines from VCF header
#zless $INPUT_VCF | grep -vFf 02_infos/unplaced_contigs_lines.txt > $OUTPUT_VCF


#sed "/^#/d" 