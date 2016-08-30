#!/bin/bash

REF=$1
REGEX_FILTER=$2

REF_MAP=$3
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  
  # setup reference
  bash A02_setup_reference.sh $REF $REF_MAP
  # log
  echo -e "$CURR_WD/$FN:bash A02_setup_reference.sh $REF $REF_MAP" >>$CURR_WD/$FN/logs/log.txt
  
done
