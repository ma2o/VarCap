#!/bin/bash

REGEX=$1
MAX_READS=8000000
for file in $( ls -d */ | grep -E "$REGEX" ); do 
  BS=$( cat $file/info.txt | grep -e 'BothS' | cut -d' ' -f2 )
  echo "$file"
  echo "FILTER:$BS"
  if [ $BS -gt $MAX_READS ]; then
    BS=$MAX_READS
  fi
  echo "READS:$BS"
  sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$BS"'/' $file/variant.config
done

