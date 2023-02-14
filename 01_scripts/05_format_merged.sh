#!/bin/bash

# Format merged output to prepare for merging with SVs from long reads
# As of 20230214, can be executed only in a jupyterhub on valeria

# srun -c 1 -p ibis_small -J 05_format_merged -o log/05_format_merged_%j.log /bin/sh 01_scripts/05_format_merged.sh &
 
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


# LOAD REQUIRED MODULES
module load R/4.1
module load bcftools 

# 1. Format header
## Extract lines for fields other than INFO, FORMAT and bcftools commands
bcftools view -h $MERGED_VCF | grep -E '=END,|SVLEN|SVTYPE|=SUPP,|=SUPP_VEC,|=GT,|ALT=<|FILTER=<|##file' | grep -E -v 'contig|bcftools|cmd' > $MERGED_UNION_DIR/VCF_lines.txt
## Contigs lines, excluding regions/contigs in $REGIONS_EX
bcftools view -h $MERGED_VCF | grep 'contig' | grep -vFf $REGIONS_EX > $MERGED_UNION_DIR/VCF_chrs.txt

## Cat these together
cat $MERGED_UNION_DIR/VCF_lines.txt $MERGED_UNION_DIR/VCF_chrs.txt > $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".header

# 2. Format VCF : add explicit alternate sequence when possible
Rscript 01_scripts/format_merged.R $MERGED_VCF $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".vcf.tmp $GENOME

# 3. Add header
cat $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".header $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".vcf.tmp > $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)"_formatted.vcf

# 4. Sort 
bcftools sort $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)"_formatted.vcf > $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".sorted.vcf

# Clean up 
rm $MERGED_UNION_DIR/VCF_lines.txt
rm $MERGED_UNION_DIR/VCF_chrs.txt
rm $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".header
rm $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)".vcf.tmp
rm $MERGED_UNION_DIR/"$(basename -s .vcf $MERGED_VCF)"_formatted.vcf