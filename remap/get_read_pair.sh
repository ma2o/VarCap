#!/bin/bash

QU1=$1
TA2=$2

zcat $QU1 | awk 'NR%4==1' | while read -r line; do
  echo -e "$line"
  PAIR2=$( zcat "$TA2" | grep -A3 -e "$line" )
  HEAD2=$( echo -e "$PAIR2" | head -n 1 )
  echo $HEAD2
  HEAD1=${HEAD2/\/2/\/1}
  echo -e "$HEAD1" 
done

