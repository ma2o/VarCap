#!/bin/bash
#$-q all.q@cube[ab]*

INFILE=$1

while read line
  do
  # extract Chr,Pos
  if [[ $line != \#* ]];then
    echo $line | awk '{ print $1" "$2" "($2+1); }'
  fi
done < $INFILE


