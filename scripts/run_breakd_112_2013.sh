#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCbd
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=breakd-%A_%a.out
#SBATCH --error=breakd-%A_%a.err

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

cd $PATH_BREAKD_DATA

# first: sam to bam conversion and indexing
# java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/SortSam.jar SO=coordinate INPUT=$SAM OUTPUT=$OUTDIR/$SAM_BASE.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# replace/add read group header
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/picard.jar AddOrReplaceReadGroups I=$SAM O=$OUTDIR/$SAM_BASE.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
# mv $OUTDIR/$SAM_BASE.rgh.bam $OUTDIR/$SAM_BASE.bam

# 1st: prepare config file
PATH=${PATH_SAMTOOLS}:$PATH
perl -I $PATH_BAM2CONF -I $PATH_GD_GRAPH_HISTOGRAM $PATH_BAM2CONF/bam2cfg_2014.pl $OUTDIR/$SAM_BASE.bam >$OUTDIR/$SAM_BASE.conf
# substitute the min_insert_size to be not shorter than the readlength -> lower:READLENGTH
sed -i "s/lower\:[0-9]*\.*[0-9]*/lower\:$MIN_INSERT_SIZE/" $OUTDIR/$SAM_BASE.conf

# 2d: run breakdancer max on config file ( parameters optimized for >100bp illumina reads)
$PATH_BREAKD/breakdancer-max -s 15 -q 35 -r 8 $PATH_BREAKD_DATA/$OUTDIR/$SAM_BASE.conf >$OUTDIR/$SAM_BASE.ctx

# postfilter output file due to cutoff values (last column of ctx file), remove is value does not exist, NA or lower than 1
CUTOFF_VAL=1
while read line; do
  if [[ "$line" =~ $( echo ^#.* ) ]]
  then
    # either print header
    echo "$line"
  else
    # or evaluate last column
    SCORE_VAL=$( echo "$line" | cut -f 12 );
    if [ $SCORE_VAL != NA ]
    then
      # echo "score_value:"$SCORE_VAL
      # if score is not NA or above 1 print out result
      if [ $( echo $SCORE_VAL">"$CUTOFF_VAL | bc -l ) -eq 1 ]
      then
        echo "$line"
      fi
    fi
  fi
done <$OUTDIR/$SAM_BASE.ctx >$OUTDIR/$SAM_BASE.filter
mv $OUTDIR/$SAM_BASE.filter $OUTDIR/$SAM_BASE.ctx

