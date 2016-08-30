#!/bin/bash

REGEX_FILTER=$1
READ_LEN=125
IS_SIZE=250
SUB_SIZE=-1
READ1_TRIM=0
READ1_TRIMF=9
REF=/proj/genomes/ctr_wga_harris/Ctr_D_UW3.fasta
CONTIG=NC_000117.fasta

for file in $( ls -d */ | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  sed -i 's#^READLENGTH=.*$#READLENGTH='"${READ_LEN}"'#' variant.config
  # sed -i 's#^INSERT_SIZE=.*$#INSERT_SIZE='"${IS_SIZE}"'#' variant.config
  sed -i 's#^SUBSAMPLE_SIZE_ALT=.*$#SUBSAMPLE_SIZE_ALT='"${SUB_SIZE}"'#' variant.config
  # sed -i 's#^READS1_TRIM=.*$#READS1_TRIM='"${READ1_TRIM}"'#' variant.config
  # sed -i 's#^READS1_TRIMF=.*$#READS1_TRIMF='"${READ1_TRIMF}"'#' variant.config
  # sed -i 's#^BAM_NAME_RAW#BAM_NAME_BASE#' $file/variant.config
  # sed -i 's#^REF_FA_ORIGINAL=.*#REF_FA_ORIGINAL='"${CONTIG}"'#' $file/variant.config
done

