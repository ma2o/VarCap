#!/bin/bash

REGEX_FILTER=$1
REF=$2
REF_MAP=$3
CURR_WD=$( pwd | sed 's/\/$//')


### 1. Set reference
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index|pdf' | grep -E "${REGEX_FILTER}" | grep -v "_99" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  # compress files
  gzip filter/*.fastq 2>/dev/null || echo "Files already zipped."
  # concatenate fwd and rev reads
  echo -e "Concatenate reads."
  if [ ! -f "filter/$( ls filter/ | grep -e "orig.gz" | head -n 1 )" ]; then
    ls filter/* | grep -e "gz$" | while read line; do
    FFN=$( echo $line | sed 's/\/$//' ); 
    echo $FFN; 
    mv $FFN $FFN.orig.gz;
    done
    zcat filter/*orig.gz >filter/${FN}_alt_1.fastq
    gzip filter/${FN}_alt_1.fastq
  fi
  mkdir -p reference
  # setup reference
  sbatch A02_search_ref.sh $REF $REF_MAP
  # log
  echo -e "$CURR_WD/$FN:bash A02_setup_reference.sh $REF $REF_MAP" >>$CURR_WD/$FN/logs/log.txt
  cd $CURR_WD
done



