#!/bin/bash
#$-q all.q@cube[ab]*
# resamples subsets of reads using subsample.py script
SUBSAMPLE_SIZE_ALT=$1
BAM_NAME_BASE=$2
REPEATS=$3
# set absolute paths with names of reference and variant reads
VAR1=$4
VAR2=$5
PATH_SCRIPTS=$6
IT=$7

VAR1_NAME=$(basename $VAR1 | sed 's/\.f.*q*$//')
VAR2_NAME=$(basename $VAR2 | sed 's/\.f.*q*$//')

READCN=$(echo $SUBSAMPLE_SIZE_ALT | sed 's/...$/k/')

SUBSETS=subsets

mkdir $SUBSETS
  
#generate unzipped files
# PATH_TEMP=/temp/varcap/subsample/${BAM_NAME_BASE}
# mkdir -p $PATH_TEMP
echo "VAR1:$VAR1"
if [[ "$VAR1" == *gz ]]; then
    echo "VARCAP: Unzipping files."
    zcat $VAR1 >$SUBSETS/$VAR1_NAME.unzip.fastq
    zcat $VAR2 >$SUBSETS/$VAR2_NAME.unzip.fastq
    VAR1=$SUBSETS/$VAR1_NAME.unzip.fastq
    VAR2=$SUBSETS/$VAR2_NAME.unzip.fastq
fi
  
# for ((  i = 1 ; i <= $REPEATS;  i++  )); do
    echo "variant subsample run: $IT"
    python $PATH_SCRIPTS/filter/subsample.py -pe -n $SUBSAMPLE_SIZE_ALT -in $VAR1 $VAR2 -out ${BAM_NAME_BASE}_r${IT}_1.fq ${BAM_NAME_BASE}_r${IT}_2.fq
# done

rm $SUBSETS/$VAR1_NAME.unzip.fastq
rm $SUBSETS/$VAR2_NAME.unzip.fastq
