#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCpi25
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=pindel25-%A_%a.out
#SBATCH --error=pindel25-%A_%a.err

. /etc/profile

BAM=$1
BAM_BASE=$(basename $BAM | sed 's/\..am$//')
OUTPUT=$2
PATH_PINDEL=$3
PATH_SAM2PINDEL=$4
PATH_REF=$5
PATH_SAMTOOLS=$6
PATH_PINDEL_DATA=$7
INSERT_SIZE=$8
REF_ID=$(sed -n '1p' $PATH_REF | sed 's/ .*$//' | sed 's/^>//')
# echo "REF_ID: "$REF_ID

cd $PATH_PINDEL_DATA
# first run sam2pindel to create pindel input file 
# bam input: ./samtools view input.bam | ./sam2pindel - Output4Pindel.txt 300 tumour 0
$PATH_SAMTOOLS/samtools view $BAM | $PATH_PINDEL/sam2pindel - $OUTPUT/$BAM_BASE.pindel.txt $INSERT_SIZE E25_varcap 0 Illumina-PairEnd

# run pindel with pindel input file (after sam2pindel) (default:-n/--min_NT_size 50 -v/--min_inversion_size 50)
$PATH_PINDEL/pindel -T 2 -f $PATH_REF -p $OUTPUT/$BAM_BASE.pindel.txt -c ALL -o $OUTPUT/$BAM_BASE

# run pindel2vcf for D and SI files
$PATH_PINDEL/pindel2vcf -G -p $OUTPUT/${BAM_BASE}_D -r $PATH_REF -R $REF_ID -d 20130121 -v $OUTPUT/${BAM_BASE}_D.vcf
$PATH_PINDEL/pindel2vcf -G -p $OUTPUT/${BAM_BASE}_SI -r $PATH_REF -R $REF_ID -d 20130121 -v $OUTPUT/${BAM_BASE}_SI.vcf
$PATH_PINDEL/pindel2vcf -G -p $OUTPUT/${BAM_BASE}_INV -r $PATH_REF -R $REF_ID -d 20130121 -v $OUTPUT/${BAM_BASE}_INV.vcf
$PATH_PINDEL/pindel2vcf -G -p $OUTPUT/${BAM_BASE}_TD -r $PATH_REF -R $REF_ID -d 20130121 -v $OUTPUT/${BAM_BASE}_TD.vcf
$PATH_PINDEL/pindel2vcf -G -P $OUTPUT/${BAM_BASE} -r $PATH_REF -R $REF_ID -d 20130121 -v $OUTPUT/${BAM_BASE}_ALL.vcf

# run pindel with bam file (use bam file from breakdancer run)
# $PATH_PINDEL/pindel -T 2 -f $PATH_REF -i $BAM -c ALL -o $OUTPUT

