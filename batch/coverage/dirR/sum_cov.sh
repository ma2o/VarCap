#!/bin/bash

INDIR=$1

COV_TABLE="NA"

for file in $( ls $INDIR | grep -E '_10.txt$'); do
  echo "$file"
  if [ "$COV_TABLE" = 'NA' ]; then
    COV_TABLE=$( less $INDIR/$file | cut -f1,2 | head -n 5 )
    less $INDIR/$file | cut -f1,2 >cov1_temp
  else
    less $INDIR/$file | cut -f2 >cov2_temp
    paste cov1_temp cov2_temp | column -s $'\t' -t >covA_temp
    mv covA_temp cov1_temp
  fi
done
rm cov2_temp

mv cov1_temp cov_final.txt

