#!/bin/bash

INDIR=$( pwd | sed 's/\/$//' )
REGEX=$1
SOURCE_DIR=$2

for file in $( ls | grep -E "$REGEX" ); do
  echo $file
  TLD=$( pwd )
  # mkdir -p $file/raw
  cd $file/raw
  FILENM=$( ls $SOURCE_DIR | grep -E "$file" )
  i=1
  for reads in ${FILENM}; do
    READS_OUT=$( echo $reads | sed 's#.fastq.gz$#_'"${i}"'.fastq.gz#' )
    echo $READS_OUT
    ln -s $SOURCE_DIR/$reads $READS_OUT
    let i++
  done
  cd $TLD
done

