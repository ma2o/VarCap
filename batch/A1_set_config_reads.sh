#!/bin/bash

# INDIR=$1
REGEX_FILTER=$1
READ_LEN=120
IS_SIZE=250
READ1_TRIM=0
READ1_TRIMF=9
REF=/proj/evochlamy/ref_upd/NC_005861_PAC_03_mito.fasta

for file in $( ls -d */ | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  echo "$FN"
  # create projects
  # cd $file
  # PROJ_NAME=$( grep -e 'PROJ_NAME=' ${file}/variant.config )
  # echo $PROJ_NAME
  sed -i 's#^BAM_NAME_RAW=.*$#BAM_NAME_RAW='"${FN}"'_all#' ${file}/variant.config
done

