#!/bin/bash

# SLURM
#SBATCH --job-name=VCsetref
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --output=log/setref-%A_%a.out
#SBATCH --error=log/setref-%A_%a.err

# Script creates necessary files and folders for the main project run.

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG


REF_ORI=$1
REF_MAP=$2
REFN=NONE
REFN_FULL=NONE

### 0. If reference is supplied, update variant config
  if [[ ! -z "$REF_ORI" ]]; then
    sed -i 's#^REF_FA_ORIGINAL=.*#REF_FA_ORIGINAL='"${REF_ORI}"'#' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
  elif [[ ! -z "$REF_MAP" ]]; then
    sed -i 's#^REF_MAPPING=.*#REF_MAPPING='"${REF_MAP}"'#' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
  fi
echo "Reference in variant.config added."
echo -e "Reference in variant.config added: ${REF_ORI} ${REF_MAP}" >>$PATH_LOGS/log.txt


### 1. Check if reference is supplied, else take from variant.config
if [[ -z "$REF_ORI" ]]; then
  REF_ORI=$REF_FA_ORIGINAL
fi

REFN=$( basename $REF_ORI | sed 's/\.f.*a$/\.fasta/')
DIRNAME=$( dirname $REF_FA )
REFN_FULL="${DIRNAME}/${REFN}"

if [[ -z "$REF_MAP" ]]; then
  REF_MAP=$REF_MAPPING
fi

### check if reference exists, else quit
if [ -e "${REFN_FULL}" ] && [ -e "${REFN_FULL}.list" ]; then
  echo "VARCAP:Reference exists, exiting: ${REFN_FULL} and ${REFN_FULL}.list"
  echo "VARCAP:Reference exists, exiting: ${REFN_FULL} and ${REFN_FULL}.list" >>$PATH_LOGS/log.txt
else
  echo "VARCAP:Processing reference: ${REFN_FULL}"
  ### 2. prepare reference list
  mkdir -p $REF_FA_MOD
  ln -s "$REF_ORI" "$REF_FA"
  # cp $REF_ORI $REF_FA_MOD
  # create lookup list for char exchange | to _ then modify fasta descriptions
  grep -e '>' $REFN_FULL | cut -d ' ' -f 1 | cut -d '>' -f 2 >${REFN_FULL}.list1
  # extract identifier from NCBI old fasta headers
  while read line; do
    if [[ $line == gi\|* ]]; then
      echo $line | cut -d '|' -f 4  >>${REFN_FULL}.list2 # sed 's/\./_/g'
    else
      echo $line  >>${REFN_FULL}.list2 # sed 's/\./_/g'
    fi
  done <${REFN_FULL}.list1
  paste -d"\t" ${REFN_FULL}.list2 ${REFN_FULL}.list1 >${REFN_FULL}.list
  rm ${REFN_FULL}.list1 ${REFN_FULL}.list2
fi

### check for mapping reference
if [[ -d $REF_MAP ]]; then
  for file in $( ls $REF_MAP | grep -E 'fna$|fasta$|fa$' ); do
    MAPREFPATHFULL=$REF_MAP/$file
    MAPREFLOCAL=$( echo "$REF_FA_MOD/$file" | sed 's/\.f.*a$/\.fasta/' )
    if [ -e "$MAPREFLOCAL" ] && [ -e "${MAPREFLOCAL}.list" ]; then
      echo "VARCAP:Mapping Reference exists, exiting: ${MAPREFLOCAL} and ${MAPREFLOCAL}.list"
    else
      mkdir -p $REF_FA_MOD
      ln -s "$MAPREFPATHFULL" "$MAPREFLOCAL"
      
      # create lookup list for char exchange | to _ then modify fasta descriptions
      grep -e '>' $MAPREFLOCAL | cut -d ' ' -f 1 | cut -d '>' -f 2 >${MAPREFLOCAL}.list1
      # extract identifier from NCBI old fasta headers
      while read line; do
        if [[ $line == gi\|* ]]; then
          echo $line | cut -d '|' -f 4  >>${MAPREFLOCAL}.list2 # sed 's/\./_/g'
        else
          echo $line  >>${MAPREFLOCAL}.list2 # sed 's/\./_/g'
        fi
      done <${MAPREFLOCAL}.list1
      paste -d"\t" ${MAPREFLOCAL}.list2 ${MAPREFLOCAL}.list1 >${MAPREFLOCAL}.list
      rm ${MAPREFLOCAL}.list1 ${MAPREFLOCAL}.list2
    fi
  done
fi


echo "Reference preperation finished."
echo -e "Reference preperation finished: ${REFN_FULL}" >>$PATH_LOGS/log.txt
