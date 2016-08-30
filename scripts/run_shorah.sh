#!/bin/bash
#$-q all.q@cube[ab]*
REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_SHORAH=$4
PATH_SHORAH=/scratch/zojer/programs/shorah-0.8
PATH_SAMTOOLS=$5

cd $OUTDIR
# sort bam file
$PATH_SAMTOOLS/samtools sort $BAM ${BAM_BASE}_sorted

#run shorah
$PATH_SHORAH/shorah.py --bam=${BAM_BASE}_sorted.bam --fasta=$REF


