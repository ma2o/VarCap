#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

for file in $( ls $INDIR -d */ | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file
  bash 400_collect_raw_variant_files.sh
  bash 401_calls2vcf_raw.sh
  bash 402_repeats2vcf.sh
  bash 403_snp_regions2vcf.sh
  bash 404_filter2vcf.sh 2
  
done


