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

# run muTect 114: needs I:tumor and I:normal samples
java -Xmx4g -jar muTect-1.1.4.jar -T MuTect -nt 2 -nct 1 -I:tumor /scratch/zojer/projects/varcap_data/09_135Av07_2M_Proto_am_U25/mapper/bwa/09_135Av07_2M_Proto_am_U25_400x16_bwa_3_v1.bam -R /scratch/zojer/projects/genomes/NC_005861/NC_005861.fasta

#run lofreq
$PATH_LOFREQ2/lofreq call -E -f $REF -o ${BAM_BASE}_lofreq2.vcf $BAM

# filter
$PATH_LOFREQ2/lofreq filter --strandbias holm-bonf --min-cov 8 -i ${BAM_BASE}_lofreq2.vcf -o ${BAM_BASE}_lofreq2_filter.vcf

