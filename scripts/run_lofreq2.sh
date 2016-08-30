#!/bin/bash

# SLURM
#SBATCH --job-name=wVCClf2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=lofreq2-%A_%a.out
#SBATCH --error=lofreq2-%A_%a.err

. /etc/profile

REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_LOFREQ2=$4
PATH_SAMTOOLS=$5

# PATH_LOFREQ2=/home/user/zojer/programs/lofreq_star-2.0.0-beta/lofreq

cd $OUTDIR
# sort bam file
# $PATH_SAMTOOLS/samtools sort $BAM $BAM_BASE.sorted

#run lofreq
$PATH_LOFREQ2/lofreq call -E -f $REF -o ${BAM_BASE}_lofreq2.vcf $BAM

# filter
$PATH_LOFREQ2/lofreq filter --strandbias holm-bonf --min-cov 8 -i ${BAM_BASE}_lofreq2.vcf -o ${BAM_BASE}_lofreq2_filter.vcf

