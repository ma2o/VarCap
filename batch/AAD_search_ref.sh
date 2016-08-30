#!/bin/bash

REFPATH1=$1
MAPPATH=$2
REGEX_FILTER=$3

CURR_WD=$( pwd | sed 's/\/$//')

for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index.*' | grep -E "${REGEX_FILTER}" | grep -Ev 'ERR108272|ERR108273' ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN

  # setup reference
  # bash $CURR_WD/$FN/A02_search_ref.sh $REFPATH1 $CURR_WD/$FN
  echo $( pwd )
  qsub -l vf=4G $CURR_WD/$FN/A02_search_ref.sh $REFPATH1 $MAPPATH
  # log
  echo -e "$CURR_WD/$FN:bash A02_setup_reference.sh $REFPATH1 $MAPPATH" >>$CURR_WD/$FN/log.txt
  cd $CURR_WD
done
