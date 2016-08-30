#!/bin/bash
#$-q all.q@cube[ab]*

BAM=$1
OUTDIR=$2
REPEATS=$3

PATH_BREAKD=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/breakdancer
PATH_PINDEL=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/pindel
PATH_DELLY=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/delly

#run i times
for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "run: $i"
  #if repeat == 1, then do not modify name and run for once, else modify ending to 1...i and run ix
  if [ $REPEATS == 1 ]
    then
      BAM_MOD=$BAM
    else
      BAM_BASE=$(basename $BAM | sed 's/.\.bam$//')
      BAM_DIR=$(dirname $BAM)
      BAM_MOD=$BAM_DIR/${BAM_BASE}${i}.bam
  fi

echo "modified: "$BAM_MOD
#run breakdancer
mkdir $PATH_BREAKD/$OUTDIR
#sh $PATH_BREAKD/run_breakd.sh $BAM_MOD $OUTDIR

#run pindel
mkdir $PATH_PINDEL/$OUTDIR
sh $PATH_PINDEL/run_pindel.sh $BAM_MOD $OUTDIR

#run delly
mkdir $PATH_DELLY/$OUTDIR
#sh $PATH_DELLY/run_delly.sh $BAM_MOD $OUTDIR

done

#check for empty files and delete them:
#check pindel dir for empty files and remove them
FILES_CHECK=$PATH_PINDEL/$OUTDIR/*
echo "Checking for empty files ..."
for f in $FILES_CHECK
do
  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
  if [ $FILESCK_SIZE == 0 ]
  then
    echo "remove "$f
    rm $f
  fi
done

#check delly dir for empty files and remove them
FILES_CHECK=$PATH_DELLY/$OUTDIR/*
echo "Checking for empty files ..."
for f in $FILES_CHECK
do
  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
  if [ $FILESCK_SIZE == 0 ]
  then
    echo "remove "$f
    rm $f
  fi
done
