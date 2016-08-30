#!/bin/bash

INFILE=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( cat $INFILE | cut -f1 | grep -v '^NAME' ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  echo $FN
  # run bwa mem
  bash C01_calling.sh
  # log
  echo -e "$CURR_WD/$FN:bash C01_calling.sh" >>$CURR_WD/$FN/logs/log.txt

done
