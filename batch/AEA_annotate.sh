#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' |grep -E "${REGEX_FILTER}"  | grep -v 'old'); do 
  echo $file;
  # create projects
  cd $CURR_WD/$file
  
  bash D01_annotate2vcf.sh
  
done


