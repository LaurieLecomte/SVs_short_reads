# SV calling pipeline from short-read sequencing data

## Pipeline Overview

1. **Prepare** required regions files in .bed and .txt : `00_prepare_regions.sh`
2. **Call SVs** : the 3 tools may be used independently in any order or at the same time.

   2.1. delly : scripts `01.1` to `01.5`
   
   2.2. manta : scripts `02.1` to `02.3`
   
   2.3. smoove : scripts `03.1` to `03.5`

3. **Merge SV calls** across callers : `04_merge_callers.sh`
4. **Format** merged output : `05_format_merged.sh`
5. **Filter** merged output : `06_filter_merged.sh` 

### Additional scripts

Other scripts targeting a specific step or operation conducted in one of the main scripts or allowing additional analyses are provided in the `01_scripts/utils` subdirectory.

* `01_scripts/utils/format_add_ALTseq.R` : adds an explicit alternate sequence (when possible) to the merged SVs. Called by the `01_scripts/utils/format_merged.R` script featured in the `05_format_merged.sh` main script.
* `01_scripts/utils/format_merged_sample_names.R` : add unique sample names to the merged VCF, also called by the `05_format_merged.sh` main script.
* `01_scripts/utils/convertInversion.py` : used for converting manta BNDs to inversions in the `02.1_manta_call` main script. It comes from the [manta GitHub repository](https://github.com/Illumina/manta/blob/75b5c38d4fcd2f6961197b28a41eb61856f2d976/src/python/libexec/convertInversion.py) and is licensed under the [GNU General Public License v3.0](https://github.com/Illumina/manta/blob/75b5c38d4fcd2f6961197b28a41eb61856f2d976/LICENSE.txt).
* `01_scripts/utils/combined_plot_by_caller.R` : used for plotting filtered short-read SVs (Supp. Fig. 2 from the paper [Investigating structural variant, indel and single nucleotide polymorphism differentiation between locally adapted Atlantic salmon populations using whole genome sequencing and a hybrid genomic polymorphism detection approach](https://www.biorxiv.org/content/10.1101/2023.09.12.557169v1))

Older scripts used for development or debugging purposes are stored in the `01_scripts/archive` folder for future reference if needed. These are not meant to be used in their current state and may be obsolete.

## Prerequisites

### Files

* A reference genome named `genome.fasta` and its index (.fai) in `03_genome`

* Bam files for all samples and their index. These can be soft-linked in the 04_bam folder for easier handling : if `$BAM_PATH` is the remote path to bam files, use `for file in $(ls -1 $BAM_PATH/*); do ln -s $file ./04_bam; done`. These should be named as `SAMPLEID.bam` (see sample ID list below). 
  * **NOTE** : a unique `RG tag` is required in bam files' entries for running `smoove` (and possibly other callers). If these read RG tags are missing or not unique, see the script [`10_replace_RG.sh`](https://github.com/enormandeau/wgs_sample_preparation/blob/master/01_scripts/10_replace_RG.sh), which adds unique RG tags *inside* bam files AND in their header.

* A bam files list in `02_infos`. This list can be generated with the following command, where `$BAM_DIR` is the path of the directory where bam files are located : `ls -1 $BAM_DIR/*.bam > 02_infos/bam_list.txt`

* A sample IDs list in `02_infos`, one ID per line. This list can be used for renaming bam files symlinks in `$BAM_DIR`, adjust `grep` command as required (warning : use carefully): `less 02_infos/ind_ALL.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bam | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam; done` and `less 02_infos/ind_ALL.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bai | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam.bai; done`

* A chromosomes list (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SV calling step. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr_list.txt`

* If some chromosomes are to be excluded from the SV calling step, such as unplaced contigs, these must be listed in `02_infos/excl_chrs.txt`, which needs to be encoded in linux format AND have a newline at the end.

* A sample IDs list (`02_infos/ind_ALL.txt`), one ID per line


### Software

#### For Manitou users
Most programs are available as modules on Manitou, which are automatically loaded in the `#LOAD REQUIRED MODULES` section in each script.
However, custom conda environments are required for running `delly` and `jasmine`, as these programs are not available on Manitou; See the [Conda environment preparation](#conda-environment-preparation) section below. 

#### For users working with other computing clusters and servers
The program versions specified in this pipeline refer to the versions available on IBIS' bioinformatics servers when this pipeline was built in 2021-2022, and are likely not available on all other servers. 
Please add a '#' at the beginning of each line in the `#LOAD REQUIRED MODULES` section in each script (or remove these lines), and follow the [Conda environment preparation](#conda-environment-preparation) to create custom conda environments with correct program versions and dependencies.
A R installation is also required.



## Detailed Walkthrough

For running each script, copy the `srun` command from the script's header to the terminal and adjust parameters (memory, partition, time limit) if necessary.  
The header also features a brief description of the script's contents. 


### Conda environment preparation

#### SV calling environment (`SVs_SR`)
From the main directory, run `conda create --name SVs_SR --file SVs_SR_env.txt`

This environment is used for calling SVs and contains the following callers:
* delly 1.1.6
* manta 1.6.0
* smoove 0.2.8 and its dependencies
* bcftools 1.13


#### SV merging environment (`jasmine_1.1.5`)
From the main directory, run `conda create --name jasmine_1.1.5 --file jasmine_1.1.5_env.txt`

This environment is used for merging SVs across callers, and contains [jasmine 1.1.5](https://github.com/mkirsche/Jasmine) and bcftools 1.13.


### Main pipeline

#### 1. Prepare region files (`00_prepare_regions.sh`)

This script prepares the bed files required for specifying the regions in which SVs must be called or must not be called. It first produces a bed file from the reference fasta in order to yield : 

* A text file of excluded chromosomes or contigs for running delly
* Bed files for each chromosome in order to parallelize manta across chromosomes
* A bed file of excluded regions, for running smoove


#### 2. Call SVs using 3 seperate tools

##### Delly (scripts 01.1 to 01.5)

Based on instructions for germline calling in high coverage genomes on [delly's GitHub](https://github.com/dellytools/delly#germline-sv-calling).

Before running each script for delly, activate the `SVs_SR` env (even if you are working on Manitou) : `conda activate SVs_SR`

* `01.1_delly_call.sh`
* `01.2_delly_merge_sites.sh`
* `01.3_delly_genotype.sh`
* `01.4_delly_merge_samples.sh`
* `01.5_delly_filter.sh`

##### Manta

* `02.1_manta_call.sh`
* `02.2_manta_merge.sh`
* `02.3_manta_filter.sh`

##### Smoove

* `03.1_smoove_call.sh`
* `03.2_smoove_merge.sh`
* `03.3_smoove_genotype.sh`
* `03.4_smoove_merge_samples.sh`
* `03.5_smoove_filter.sh`

#### 3. Merge SV calls across callers (`04_merge_callers.sh`)
Before running this script, activate the `jasmine_1.1.5` env (even if you are working on Manitou): `conda activate jasmine_1.1.5`

#### 4. Format merged output (`05_format_merged.sh`)

#### 5. Filter merged SVs (`06_filter_merged.sh`)
Keep SVs supported by at least 2/3 tools and larger than 50 bp.

