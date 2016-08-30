#!/bin/bash

REGEX_FILTER=$1
REF=$2
REF_MAP=$3
CURR_WD=$( pwd | sed 's/\/$//')


### 4. Qualty filtering
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$file
  bash A03_trimmomatic.sh
  # log
  echo -e "$CURR_WD/$FN:bash A03_trimmomatic.sh" >>$CURR_WD/$FN/log.txt
  cd $CURR_WD
done

