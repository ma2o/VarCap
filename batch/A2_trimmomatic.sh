#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

for file in $( ls $INDIR | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  bash 1101_adapter_trimmomatic.sh
  
done
