#!/bin/bash

find -maxdepth 1 -type d -name '*ERR*' | sed 's/^\.\///;s/_[0-9]*$//' | sort -u | while read -r line; do
  echo -e "$line\t"
  for file in $( ls -d */ | grep -e "$line" ); do
    echo -en "${file%%/}\t"
    echo -e "$( cat $file/reference/*summary.txt | grep -e topref | cut -f2,7 | tr "\n" "\t" )"
    # echo -en "${file%%_[0-9]*/}\t"
    # echo -e "$( cat $file/reference/*summary.txt | grep -e topref | cut -f2 )"
  done
  echo ""
done

