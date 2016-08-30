#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')


### 1. Set reference
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index|pdf' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  FNS=${FN%_*}
  # create projects
  cd $CURR_WD/$FN
  # concatenate fwd and rev reads
  echo -e "Prepare remapping."
  mv filter filter_refsearch
  cp variant.config variant.config.refsearch
  mkdir -p filter
  # mkdir -p raw
  # generate links to raw files
  cp -r $CURR_WD/../$FNS/raw .
  
  cd $CURR_WD
done



