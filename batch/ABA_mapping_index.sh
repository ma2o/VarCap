#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  
  # setup reference
  bash B01_set_mapping_index.sh
  # log
  echo -e "$CURR_WD/$FN:bash B01_set_mapping_index.sh" >>$CURR_WD/$FN/log.txt

done
