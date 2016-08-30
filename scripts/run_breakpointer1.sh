#!/bin/bash
#$-q all.q@cube[ab]*
REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_MUTECT=$4
PATH_SAMTOOLS=$5

cd $OUTDIR

# add read group header
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/AddOrReplaceReadGroups.jar I=$OUTPUT_DIR/$SAM_BASE.marked.bam O=$OUTPUT_DIR/$SAM_BASE.marked.rgroup.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
# index bam
samtools index rg_update.bam

# run breakpointer
bash BreakPointer.sh /scratch/zojer/programs/MATLAB/MATLAB_Compiler_Runtime/v83 Bp_test none /scratch/zojer/programs/muTect-1.1.4-bin/rg_update.bam none 250 6 80 200 40 50 80 15 5 10 500 100000 100 20 8 0.75

#run lofreq
$PATH_LOFREQ2/lofreq call -E -f $REF -o ${BAM_BASE}_lofreq2.vcf $BAM

# filter
$PATH_LOFREQ2/lofreq filter --strandbias holm-bonf --min-cov 8 -i ${BAM_BASE}_lofreq2.vcf -o ${BAM_BASE}_lofreq2_filter.vcf

