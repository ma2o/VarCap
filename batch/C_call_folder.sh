#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

for file in $( ls $INDIR -d */ | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  bash 31_calling.sh
  
done


