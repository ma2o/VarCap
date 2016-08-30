#!/bin/bash

REGEX_FILTER=$1

for file in $( ls -d */ | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  cd $file
  # create projects
  rm caller/cortex/ref/*
  rm caller/cortex/*
  bash 01_setup_reference.sh
  # bash 02_setup_index.sh
  
done


