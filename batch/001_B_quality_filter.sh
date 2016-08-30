#!/bin/bash

INFILE=$1
CURR_WD=$( pwd | sed 's/\/$//')

for file in $( cat $INFILE | cut -f1 | grep -v '^NAME' ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$file
  bash A03_trimmomatic.sh
  # log
  echo -e "$CURR_WD/$FN:bash A03_trimmomatic.sh" >>$CURR_WD/$FN/log.txt
  cd $CURR_WD
done

