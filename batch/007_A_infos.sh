#!/bin/bash

REGEX_FILTER=$1
CURR_WD=$( pwd | sed 's/\/$//')



### 1. get infos about filtering and mapping
echo -e "NAME\tFILTER\tTREF\tTREF_PCT\tCREF\tCREF_PCT\tFILTER_READS\tFILTER_COV\tMAPAV_COV"
for file in $( ls -d */ | grep -Ev 'reference|batch|bwa_index' | grep -E "${REGEX_FILTER}" ); do 
  FN=$( echo $file | sed 's/\/$//' )
  # create projects
  cd $CURR_WD/$FN
  # get filter infos
  NAME=$( cat info.txt | grep -e '^name' | cut -f2 )
  RL=$( cat info.txt | grep -e '^average' | head -n1 | cut -f2 )
  FILTER=$( cat info.txt | grep -e '^Both' | head -n1 | cut -d' ' -f2 )
  MAPAV=$( cat info.txt | grep -e '^map_average' | head -n1 | cut -f2 )
  # TODO cov/genome
  TREF=$( cat reference/*summary.txt | grep -e '^topref' | cut -f2 )
  TREF_LEN=$( cat reference/*summary.txt | grep -e '^topref' | cut -f4 )
  TREF_PCT=$( cat reference/*summary.txt | grep -e '^topref' | cut -f7 )
  CREF=$( cat reference/*summary.txt | grep -e '^extref' | head -n1 | cut -f2 )
  CREF_PCT=$( cat reference/*summary.txt | grep -e '^extref' | head -n1 | cut -f7 )
  if [ -z $CREF ]; then
    CREF="nd"
    CREF_PCT="0"
  fi
  # calculate theoretical mapping coverage
  FILTER_COV=$( echo -e "$FILTER\t$TREF_PCT\t$RL\t$TREF_LEN" | awk '{ reads=$1*$2*2/100;cov=reads*$3/$4; printf "%d %d", reads,cov }' )
  
  echo -e "$NAME\t$FILTER\t$TREF\t$TREF_PCT\t$CREF\t$CREF_PCT\t$FILTER_COV\t$MAPAV"
  
  # mapping stats
#   for file in $( ls -d */ | grep -E 'ERR07|ERR10|ERR11' ); do 
#	echo -en "$file\t"; 
#	cat $file/vcfs/*filter_2.vcf | grep -e SNP | grep -v 'CV1' | sort -k2,2 -un | wc -l | tr '\n' '\t';
#	cat $file/vcfs/*filter_2.vcf | grep -e SNP | grep -v 'CV1' | sort -k2,2 -un | grep -Ev '100$|99\.[0-9]*$|98\.[0-9]*$' | wc -l; 
#  done
  
  cd $CURR_WD
done
