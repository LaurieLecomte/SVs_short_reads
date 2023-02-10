#!/bin/bash

# Format SV calls by manta prior to merging across callers

# srun -c 1 -p ibis_small -J 02.4_manta_format -o log/02.4_manta_format_%j.log /bin/sh 01_scripts/02.4_manta_format.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_LIST="02_infos/chr_list.txt"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"


REGIONS_EX="02_infos/excl_chrs.txt"

MANTA_VCF="$FILT_DIR/manta/manta_PASS_PRECISE.vcf.gz"

# LOAD REQUIRED MODULES
module load bcftools
module load r/4.1.2


# 1. Construct simplified header
## Extract lines for fields other than INFO, FORMAT and bcftools commands
bcftools view -h $MANTA_VCF | grep -E -v 'INFO|FORMAT|contig|bcftools|cmd' > $FILT_DIR/manta/VCF_lines.txt
## Contigs lines, excluding regions/contigs in $REGIONS_EX
bcftools view -h $MANTA_VCF | grep 'contig' | grep -vFf $REGIONS_EX > $FILT_DIR/manta/VCF_chrs.txt

## Add line describing ALT seq INFO field
echo '##INFO=<ID=ALT_MANTA,Number=1,Type=String,Description="Alternate sequence for SV">' >> $FILT_DIR/manta/VCF_lines.txt

## Cat these together
cat $FILT_DIR/manta/VCF_lines.txt 02_infos/shared_header_lines.txt $FILT_DIR/manta/VCF_chrs.txt > $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".header

# 2. Format VCF contents 
Rscript 01_scripts/format_manta.R $MANTA_VCF $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".contents $GENOME

# 3. Concatenate header with contents
cat $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".header $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".contents | bcftools sort > $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)"_formatted.vcf

# Index 
#bgzip $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)"_formatted.vcf
#tabix -p vcf $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)"_formatted.vcf.gz

# Clean up 
rm $FILT_DIR/manta/VCF_lines.txt
rm $FILT_DIR/manta/VCF_chrs.txt
rm $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".header
rm $FILT_DIR/manta/"$(basename -s .vcf.gz $MANTA_VCF)".contents