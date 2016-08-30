#!/bin/bash

INDIR=$( pwd )
REGEX_FILTER=$1

for file in $( ls $INDIR -d */ | grep -Ev 'batch|reference|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  bash D03_ccf2mysql.sh vcfs varcap_2014
  cd ..

done


