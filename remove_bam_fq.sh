#!/bin/bash
#$-q all.q@cube[ab]*

for folder in $( ls | grep -e 'evochlamy' );
do
if [ -d $folder ];
then
  echo $folder
  cd $folder
  cp /scratch/zojer/projects/varcap_2.0/50_cleanup.sh .
  bash 50_cleanup.sh
  cd /scratch/zojer/projects/varcap_data
fi

done

