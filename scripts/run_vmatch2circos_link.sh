#!/bin/bash
#$-q all.q@cube[ab]*

INFILE=$1

CHROM_ID="NC_005861"
COUNTER="0"

while read line
  do
  # extract length and pos of pair and assign to unique id
  if [[ $line != \#* ]];then
    LENGTH1=$(echo $line | awk '{print $1;}' )
    POS1=$(echo $line | awk '{print $2;}' )
    POS1_END=$(( $POS1 + $LENGTH1 - 1 ))
    LENGTH2=$(echo $line | awk '{print $4;}' )
    POS2=$(echo $line | awk '{print $5;}' )
    POS2_END=$(( $POS2 + $LENGTH2 - 1 ))
    UNIQUE_ID="seqdup"
    echo $UNIQUE_ID$COUNTER" "$CHROM_ID" "$POS1" "$POS1_END
    echo $UNIQUE_ID$COUNTER" "$CHROM_ID" "$POS2" "$POS2_END
    COUNTER=$(($COUNTER+1))
  fi
done < $INFILE

