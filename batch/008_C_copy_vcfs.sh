#!/bin/bash

INFILE=$1
CURR_WD=$( pwd | sed 's/\/$//')
TARDIR=/proj/genomes/ctr_sanger/results_remap/vcfs

for file in $( cat $INFILE | cut -f1 | grep -v '^NAME' ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  echo $FN
  # copy vcfs to proj
  mkdir $TARDIR/$FN
  mkdir $TARDIR/$FN/vcfs
  mkdir $TARDIR/$FN/vcfs_raw
  cp *txt *config *sh $TARDIR/$FN/
  cp -r vcfs/* $TARDIR/$FN/vcfs/
  cp -r vcfs_raw/* $TARDIR/$FN/vcfs_raw/

done
