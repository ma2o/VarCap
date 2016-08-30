#!/bin/bash
#$-q all.q@cube[ab]*

. /etc/profile

SAM=$1
OUTDIR=$2
SAM_BASE=$(basename $SAM | sed 's/\..am$//')
PATH_BREAKD=$3
PATH_BAM2CONF=$4
PATH_SAMTOOLS=$5
PATH_PICARD=$6
PATH_BREAKD_DATA=$7
PATH_GD_GRAPH_HISTOGRAM=$8
MIN_INSERT_SIZE=$9
PATH_STATISTICS_DESCRIPTIVE=/scratch/zojer/programs/Statistics-Descriptive-3.0607/lib
PATH_CDF=/scratch/zojer/programs/Math-CDF-0.1

#export PERL5LIB=/home/user/zojer/programs/breakdancer-1.1_2011_02_21/breakdancer-1.1_2011_02_21/perl"${PERL5LIB:+:$PERL5LIB}"
#or
#export PERL5LIB=${PERL5LIB}:/home/user/zojer/programs/breakdancer-1.1_2011_02_21/breakdancer-1.1_2011_02_21/cpp
#or
#perl -I /home/user/zojer/programs/breakdancer-1.1_2011_02_21/breakdancer-1.1_2011_02_21/perl $PATH_BAM2CONF/bam2cfg.pl

cd $PATH_BREAKD_DATA

#first: sam to bam conversion and indexing
# java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/SortSam.jar SO=coordinate INPUT=$SAM OUTPUT=$OUTDIR/$SAM_BASE.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

#replace/add read group header
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/AddOrReplaceReadGroups.jar I=$SAM O=$OUTDIR/$SAM_BASE.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
# mv $OUTDIR/$SAM_BASE.rgh.bam $OUTDIR/$SAM_BASE.bam

#1st: prepare config file
PATH=${PATH_SAMTOOLS}:$PATH
perl -I $PATH_BAM2CONF -I $PATH_GD_GRAPH_HISTOGRAM $PATH_BAM2CONF/bam2cfg_2014.pl $OUTDIR/$SAM_BASE.bam >$OUTDIR/$SAM_BASE.conf

#2d: run breakdancer max on config file ( parameters optimized for >100bp illumina reads)
# to avoid false positives from to short insert sizes, the min insert size is the readlength -s
$PATH_BREAKD/breakdancer_max -s 15 -q 35 -r 8 -s $MIN_INSERT_SIZE $PATH_BREAKD_DATA/$OUTDIR/$SAM_BASE.conf >$OUTDIR/$SAM_BASE.ctx
# $PATH_BREAKD/breakdancer_max $PATH_BREAKD_DATA/$OUTDIR/$SAM_BASE.conf >$OUTDIR/$SAM_BASE.ctx

