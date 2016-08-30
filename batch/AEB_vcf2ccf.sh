#!/bin/bash

INDIR=$( pwd )
REGEX_FILTER=$1
LISTFILE=$2

LISTFILE=/scratch/zojer/projects/varcap_pilot_04_2/evochlamy_list.csv

for file in $( ls $INDIR -d */ | grep -Ev 'batch|reference|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  INREG=$( basename $INDIR/$file | grep -Eo '^[0-9]{5}' )
  echo $INREG
  INVAR=$( cat $LISTFILE | grep -E "$INREG" )
  echo $INVAR
  bash D02_vcf2ccf.sh $INVAR
  cd ..

done

