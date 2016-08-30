#!/bin/bash

MAP_READCOUNT=$2
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
  cd $CURR_WD
done

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  if [[ ! -z "$MAP_READCOUNT" ]]; then 
    sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$MAP_READCOUNT"'/' variant.config
  fi
  # run bwa mem
  bash B03_mapping_bwamem.sh
  # log
  echo -e "$CURR_WD/$FN:bash B03_mapping_bwamem.sh" >>$CURR_WD/$FN/log.txt
done
