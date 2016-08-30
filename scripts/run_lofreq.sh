#!/bin/bash

# SLURM
#SBATCH --job-name=wVCClf1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=lofreq1-%A_%a.out
#SBATCH --error=lofreq1-%A_%a.err

. /etc/profile

REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_LOFREQ=$4
PATH_SAMTOOLS=$5

echo "update paths"
# update paths
PATH=/apps/lofreq/0.6.1/bin:$PATH
export PYTHONPATH=/apps/lofreq/0.6.1/lib64/python2.7/site-packages/:$PYTHONPATH
echo $PATH

cd $OUTDIR
# sort bam file
# $PATH_SAMTOOLS/samtools sort $BAM $BAM_BASE.sorted

#run lofreq
$PATH_LOFREQ/lofreq_snpcaller.py -f $REF -b $BAM -o ${BAM_BASE}_lofreq

# filter
$PATH_LOFREQ/lofreq_filter.py --strandbias-holmbonf --min-cov 8 -i ${BAM_BASE}_lofreq -o ${BAM_BASE}_lofreq_filter

