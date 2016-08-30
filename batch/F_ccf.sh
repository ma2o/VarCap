#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

LISTFILE=/scratch/zojer/projects/varcap_data2/evochlamy_list.csv

for file in $( ls $INDIR | grep -E '^[0-9]' | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  INREG=$( basename $INDIR/$file | grep -Eo '^[0-9]{5}' )
  echo $INREG
  INVAR=$( less $LISTFILE | grep -E "$INREG" )
  echo $INVAR
  bash 407_vcf2ccf.sh $INVAR
  
done


