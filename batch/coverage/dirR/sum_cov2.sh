#!/bin/bash

INDIR=$1

ROWS="NA"
ARRAYZ=()

for file in $( ls $INDIR | grep -E '_10.txt$'); do
  echo "$file"
  if [ "$COV_TABLE" = 'NA' ]; then
    
    ROWS=$( less $INDIR/$file | cut -f1 | head -n 5 )
    VAL1=$( less $INDIR/$file | cut -f2 | head -n 5 )
    ${ARRAZ[0]}=$ROWS
    ${ARRAZ[1]}=$VAL1
    echo "$COV_TABLE[@]"
  else
    COV_TABLE_TEMP=$( less $INDIR/$file | cut -f2 | head -n 5 )
    
    # COV_TABLE=$COV_TABLE_NEW
    echo "$COV_TABLE"
    echo "$COV_TABLE_NEW"
  fi
  echo "$COV_TABLE"
done

