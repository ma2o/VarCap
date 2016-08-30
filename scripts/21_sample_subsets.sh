#!/bin/bash
#$-q all.q@cube[ab]*
# resamples subsets of reads using subsample.py script
SUBSAMPLE_SIZE_REF=$1
SUBSAMPLE_SIZE_ALT=$2
BAM_NAME_BASE=$3
REPEATS=$4
# set absolute paths with names of reference and variant reads
VAR1=$5
VAR2=$6
REF1=$7
REF2=$8
PATH_SCRIPTS=$9

REF1_NAME=$(basename $REF1 | sed 's/\.f.*q$//')
REF2_NAME=$(basename $REF2 | sed 's/\.f.*q$//')
VAR1_NAME=$(basename $VAR1 | sed 's/\.f.*q$//')
VAR2_NAME=$(basename $VAR2 | sed 's/\.f.*q$//')

RC_RATIO=$(((100 * ${SUBSAMPLE_SIZE_ALT}) / (${SUBSAMPLE_SIZE_REF} + ${SUBSAMPLE_SIZE_ALT}) | bc))
REF_READCN=$(echo $SUBSAMPLE_SIZE_REF | sed 's/...$/k/')
READCN=$(echo $SUBSAMPLE_SIZE_ALT | sed 's/...$/k/')

SUBSETS=subsets

mkdir $SUBSETS

if [ $SUBSAMPLE_SIZE_REF == 0 ]
then
  # create subsample of variant (assuming realdata)
  for ((  i = 1 ; i <= $REPEATS;  i++  ))
  do
    echo "variant subsample run: $i"
    python $PATH_SCRIPTS/filter/subsample.py -pe -n $SUBSAMPLE_SIZE_ALT -in $VAR1 $VAR2 -out ${BAM_NAME_BASE}_${RC_RATIO}_r${i}_1.fq ${BAM_NAME_BASE}_${RC_RATIO}_r${i}_2.fq
  done
else
  # create mix of subsampled reads for ref and variant (assuming simulated data)
  echo "creating reference subset"
  python $PATH_SCRIPTS/filter/subsample.py -pe -n $SUBSAMPLE_SIZE_REF -in $REF1 $REF2 -out ${SUBSETS}/${REF1_NAME}_${REF_READCN}_1.fq ${SUBSETS}/${REF2_NAME}_${REF_READCN}_2.fq
  # create subsample of variant, concatenate(cat) with subsample of reference and write to new read pair
  for ((  i = 1 ; i <= $REPEATS;  i++  ))
  do
    echo "subsample run: $i"
    python $PATH_SCRIPTS/filter/subsample.py -pe -n $SUBSAMPLE_SIZE_ALT -in $VAR1 $VAR2 -out ${VAR1_NAME}_$SUBSAMPLE_SIZE_ALT ${VAR2_NAME}_$SUBSAMPLE_SIZE_ALT
    cat ${SUBSETS}/${REF1_NAME}_${REF_READCN}_1.fq ${VAR1_NAME}_$SUBSAMPLE_SIZE_ALT >${BAM_NAME_BASE}_${RC_RATIO}_r${i}_1.fq
    cat ${SUBSETS}/${REF2_NAME}_${REF_READCN}_2.fq ${VAR2_NAME}_$SUBSAMPLE_SIZE_ALT >${BAM_NAME_BASE}_${RC_RATIO}_r${i}_2.fq
    rm ${VAR1_NAME}_$SUBSAMPLE_SIZE_ALT ${VAR2_NAME}_$SUBSAMPLE_SIZE_ALT
  done
fi

