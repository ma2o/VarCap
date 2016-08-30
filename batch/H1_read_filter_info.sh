#!/bin/bash

REGEX=$1

for file in $( ls -d */ | grep -E "$REGEX" ); do 
  echo -n "$file ";
  cat $file/filter/run_trimmomatic.sh.e* | grep -E "Input Read Pairs" | sed 's/ Forward Only.*//' | sed 's/ *//g' | sed 's/Both/\tBoth/g';
  cat $file/filter/run_trimmomatic.sh.e* | grep -E "Input Read Pairs" | sed 's/ Forward Only.*//' | sed 's/ *//g' | sed 's/\:/ /g' | sed 's/Both/\nBoth/g' | sed 's/(/ (/' >>$file/info.txt
done

