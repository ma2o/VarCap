#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCdl72
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --nice=1000
#SBATCH --output=delly72-%A_%a.out
#SBATCH --error=delly72-%A_%a.err

. /etc/profile

BAM=$1
BAM_NAME=$(basename $BAM | sed 's/\..am$//')
OUTDIR=$2
REF_FA=$3
REF_FA_BASE=$( basename $REF_FA )
PATH_DELLY=$4
PATH_DELLY_DATA=$5
INSERT_SIZE=$6
REF_FA_MOD=$7
INSERT_SIZE_CUTOFF=$( echo $(( $INSERT_SIZE + $INSERT_SIZE / 3 )) )
REF=$REF_FA_MOD/$REF_FA_BASE
#
cd $PATH_DELLY_DATA/$OUTDIR
# delly version 72 parallel
$PATH_DELLY/delly_v0.7.2_parallel_linux_x86_64bit -t DEL -i 2000 -o $BAM_NAME.del.vcf -g $REF $BAM
$PATH_DELLY/delly_v0.7.2_parallel_linux_x86_64bit -t INS -i 2000 -o $BAM_NAME.ins.vcf -g $REF $BAM
$PATH_DELLY/delly_v0.7.2_parallel_linux_x86_64bit -t DUP -o $BAM_NAME.dup.vcf -g $REF $BAM
$PATH_DELLY/delly_v0.7.2_parallel_linux_x86_64bit -t TRA -o $BAM_NAME.tra.vcf -g $REF $BAM
$PATH_DELLY/delly_v0.7.2_parallel_linux_x86_64bit -t INV -o $BAM_NAME.inv.vcf -g $REF $BAM


