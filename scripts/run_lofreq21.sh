#!/bin/bash

# SLURM
#SBATCH --job-name=wVCClf21
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=lofreq21-%A_%a.out
#SBATCH --error=lofreq21-%A_%a.err

. /etc/profile

REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_LOFREQ21=$4
PATH_SAMTOOLS=$5

cd $OUTDIR

# add indel qualities
$PATH_LOFREQ21/lofreq indelqual --dindel -f $REF -o ${BAM_BASE}_lofreq21.bam $BAM

#run lofreq
$PATH_LOFREQ21/lofreq call --call-indels --min-cov 6 -f $REF -o ${BAM_BASE}_lofreq21.vcf ${BAM_BASE}_lofreq21.bam

# filter
# $PATH_LOFREQ21/lofreq filter --sb-mtc fdr --sb-incl-indels --cov-min 8 -i ${BAM_BASE}_lofreq21.vcf -o ${BAM_BASE}_lofreq21_filter.vcf

