#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCfb
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=freebayes-%A_%a.out
#SBATCH --error=freebayes-%A_%a.err

source /etc/profile

REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_FREEBAYES=$4

cd $OUTDIR

# run freebayes
echo "run freebayes pooled"
$PATH_FREEBAYES/freebayes -f $REF -F 0.01 -C 5 --pooled-continuous $BAM >${BAM_BASE}_pool.vcf
