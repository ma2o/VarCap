#!/bin/bash
#$-q all.q@cube[ab]*

#changes read name endings to new suffix
IN1=$1
OUT1=$(basename $IN1)
NEW_SUFFIX=$2

#ABC="thisvar is the rest var"
#echo ${ABC%% *}

# loop through file using awk -F field seperator, NR line number,
awk -F' ' '{
if ($1 ~ /^@/ && 0 == (NR+1) % 2) 
  print $1"'"$NEW_SUFFIX"'";
else
  print $0;
}' $IN1 >readnm_corr_$OUT1

