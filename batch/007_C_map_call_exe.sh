#!/bin/bash

# report all bam files that are smaller than Mb or Gb, which indicate a faulty mapping
for file in $( ls -d */ | grep -e ERR1 ); do echo -en "$( echo -e "$file" | sed 's/\/$//')\t"; echo -e "$( ll -h $file/mapper/bwa/*1.bam | cut -d' ' -f5 )"; done | grep -Ev "M$|G$" >not_mapped

# corect files and entries for mapping to work
cat not_mapped.txt | while read -r line; do 
  FN=$( echo -e "$file" | cut -f1 );
  echo -en "$FN\t";
  rm "$FN/filter/*fastq\*";
  sed -i 's/SUBSAMPLE_SIZE_ALT=0/SUBSAMPLE_SIZE_ALT=-1/' $FN/variant.config;
  
done
