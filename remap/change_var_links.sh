#!/bin/bash

for file in $( find . -maxdepth 1 -type d -name '*_99' | sed -e 's|^.*\/||' ); do 
  echo $file; 
  sed -i 's#^PATH_ALT_READS1=.*#PATH_ALT_READS1=${PATH_PROJECTS_DATA}/${PROJ_NAME}/filter/${PROJ_NAME}_alt_1.fastq.gz#;s#^PATH_ALT_READS2=.*#PATH_ALT_READS2=${PATH_PROJECTS_DATA}/${PROJ_NAME}/filter/${PROJ_NAME}_alt_2.fastq.gz#' $file/variant.config; 
  
done

