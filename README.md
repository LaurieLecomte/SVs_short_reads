# SV calling pipeline from short read sequencing data

## TO DO
* Load required modules outside of scripts and correct versions for both manitou and valeria 

## Pipeline Overview

1. Prepare required regions files in .bed and .txt : `00_prepare_regions.sh`
2. **Call SVs** : the 3 tools may be used independently in any order or at the same time.
* 2.1. delly : scripts `01.1` to `01.5`
* 2.2. manta : scripts `02.1` to `02.3`
* 2.2. smoove : scripts `03.1` to `03.5`
3. **Merge SV calls** across callers : `04_merge_callers.sh`
4. Format merged output : `05_format_merged.sh`
5. **Filter** merged output : `06_filter_merged.sh` 


## Prerequisites

### Files

* A reference genome named `genome.fasta` and its index (.fai) in `03_genome`
* Bam files for all samples and their index. These can be soft-linked in the 04_bam folder for easier handling : if `$BAM_PATH` is the remote path to bam files, use `for file in $(ls -1 $BAM_PATH/*); do ln -s $file ./04_bam; done`. These should be named as `SAMPLEID.bam` (see sample ID list below).
* A bam files list in `02_infos`. This list can be generated with the following command, where `$BAM_DIR` is the path of the directory where bam files are located : `ls -1 $BAM_DIR/*.bam > 02_infos/bam_list.txt`
* A sample IDs list in `02_infos`, one ID per line. This list can be used for renaming bam files symlinks in `$BAM_DIR`, adjust `grep` command as required (warning : use carefully): `less 02_infos/ind_ALL.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bam | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam; done` and `less 02_infos/ind_ALL.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bai | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam.bai; done`
* A chromosomes list (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SV calling step. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr_list.txt`
* If some chromosomes are to be excluded from the SV calling step, such as unplaced contigs, these must be listed in `02_infos/excl_chrs.txt`, which needs to be encoded in linux format AND have a newline at the end.
* A sample IDs list (`02_infos/ind_ALL.txt`), one ID per line


### Software

#### For Manitou and Valeria users



## Detailed Walkthrough

### `00_prepare_regions.sh`

This script prepares the bed files required for specifying the regions in which SVs must be called or must not be called. It first produces a bed file from the reference fasta in order to yield : 

* A text file of excluded chromosomes or contigs for running delly
* Bed files for each chromosome in order to parallelize manta across chromosomes
* A bed file of excluded regions, for running smoove


### Delly

Based on instructions for germline calling in high coverage genomes on [delly's GitHub](https://github.com/dellytools/delly#germline-sv-calling).

#### `01.1_delly_call.sh`

Call SVs in all samples in regions other than the ones in `02_infos/excl_chrs.txt`.

* On manitou : `parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -c 1 --mem=20G -p medium --time=7-00:00 -J 01.1_delly_call_{} -o log/01.1_delly_call_{}_%j.log 01_scripts/01.1_delly_call.sh {} &`
* On valeria : `parallel -a 02_infos/ind_ALL.txt -k -j 10 srun -c 1 --mem=20G -p ibis_medium --time=7-00:00 -J 01.1_delly_call_{} -o log/01.1_delly_call_{}_%j.log 01_scripts/01.1_delly_call.sh {} &`

#### `01.2_delly_merge_sites.sh`
#### `01.3_delly_genotype.sh`
#### `01.4_delly_merge_samples.sh`
#### `01.5_delly_filter.sh`

### Manta

#### `02.1_manta_call.sh`

**NOTE** : The script `convertInversion.py` comes from the manta GitHub repository (https://github.com/Illumina/manta/blob/75b5c38d4fcd2f6961197b28a41eb61856f2d976/src/python/libexec/convertInversion.py) and is licensed under the GNU General Public License v3.0 (https://github.com/Illumina/manta/blob/75b5c38d4fcd2f6961197b28a41eb61856f2d976/LICENSE.txt)


#### `02.2_manta_merge.sh`
#### `02.3_manta_filter.sh`

### Smoove

#### `03.1_smoove_call.sh`
#### `03.2_smoove_merge.sh`
#### `03.3_smoove_genotype.sh`
#### `03.4_smoove_merge_samples.sh`
#### `03.5_smoove_filter.sh`

### `04_merge_callers.sh`

### `05_format_merged.sh`

### `06_filter_merged.sh`


