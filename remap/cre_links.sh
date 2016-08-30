#!/bin/bash

for file in $( ls -d */ | grep -e "ERX321689_ERR348853" ); do 
  FN=${file%_*/}; 
  echo $FN; 
  cd $file; 
  mkdir raw; 
  cd raw;
  pwd 
  ln -s /scratch/zojer/projects/sanger_ctr/ebi_hs/$FN/ERR348853_1.fastq.gz "${FN}_ERR348853_1.fastq.gz"
  ln -s /scratch/zojer/projects/sanger_ctr/ebi_hs/$FN/ERR348853_2.fastq.gz "${FN}_ERR348853_2.fastq.gz"
  cd ../..; 
done


