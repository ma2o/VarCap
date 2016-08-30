#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$file
  bash A03_trimmomatic.sh
  # log
  echo -e "$CURR_WD/$FN:bash A03_trimmomatic.sh" >>$CURR_WD/$FN/log.txt
done
