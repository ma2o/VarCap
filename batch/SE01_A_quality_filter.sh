#!/bin/bash

REF=$1
REGEX_FILTER=$2
REF_MAP=$3
CURR_WD=$( pwd | sed 's/\/$//')

# 1. Check if reference was set, else report and exit
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do
  cd $CURR_WD/$file
  REF_SET=$( cat variant.config | grep -e '^REF_FA_ORIGINAL=' )
  if [[ "$REF_SET" =~ '/path/to/genome.fasta' && ! -f $REF ]]; then
    echo "Reference not set, please supply reference in fasta format"
    echo "Usage: bash 001_A_quality_filter.sh <reference>"
    exit 1
  fi
  cd $CURR_WD
done

### 1. Set reference
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  
  # setup reference
  bash A02_setup_reference.sh $REF $REF_MAP
  # log
  echo -e "$CURR_WD/$FN:bash A02_setup_reference.sh $REF $REF_MAP" >>$CURR_WD/$FN/logs/log.txt
  cd $CURR_WD
done


 

### 3. Run quality control on raw data
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

### 4. Qualty filtering
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$file
  bash A03_prinseq-lite.sh
  # log
  echo -e "$CURR_WD/$FN:bash A03_prinseq-lite.sh" >>$CURR_WD/$FN/log.txt
  cd $CURR_WD
done

