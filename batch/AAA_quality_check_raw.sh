#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  echo "$FN"
  # run quality checks for raw data
  INFILE_RAW=($( find  $FN/raw/* | grep -E 'fq$|fastq$|gz$|bam$'))
  for line in "${INFILE_RAW[@]}"; do
    # echo -e "$CURR_WD/$line"
    cd $FN
    qsub -l vf=4G A01_quality_check_raw.sh $CURR_WD/$line
    # log file
    echo -e "qsub $FN/A01_quality_check_raw.sh $CURR_WD/$line" >>$CURR_WD/$FN/log.txt
    cd $CURR_WD
  done
  
done

