#!/bin/bash

INDIR=$1
REGEX_FILTER=$2

for file in $( ls $INDIR | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  # create projects
  cd $INDIR/$file/mapper
  rm subsample/subsets/*
  rm subsample/*
  
done


