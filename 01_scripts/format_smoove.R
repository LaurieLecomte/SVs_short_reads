# Format smoove VCF 


# 0. Access to files provided in command line arguments -------------------
argv <- commandArgs(T)
VCF <- argv[1]            # original, unformatted vcf, can be gzipped
formatted_VCF <- argv[2]  # output vcf, will be overwritten
GENOME <- argv[3]         # reference fasta

# 1. Source required functions --------------------------------------------
source('01_scripts/add_ALT_format_smoove.R')

# 2. Process delly VCF ----------------------------------------------------
add_ALT(input_vcf = VCF, 
        output_vcf = formatted_VCF,
        refgenome = GENOME)


