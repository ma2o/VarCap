#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

for file in $( ls $INDIR | grep -E '^[0-9]' | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  bash 408_ccf2mysql.sh vcfs varcap_2014
  
done


