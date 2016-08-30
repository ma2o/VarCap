#!/bin/bash

for file in $( find -maxdepth 1 -type d -name '*ERR*' ); do
  echo -e ${file#./}
  # gzip ${file#./}/filter/*alt_1.fastq --suffix=_joined.gz
  for file2 in $( find ${file#./}/filter -name '*orig.gz' ); do
    mv $file2 ${file2%.orig.gz}
  done
done

