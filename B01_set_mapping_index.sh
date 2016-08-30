#!/bin/bash

# Script creates necessary files and folders for the main project run.

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

mkdir -p $PATH_BWA_DATA

REFNAME=$( echo $REF_FA_NAME | sed 's/\.f.*a$//')
# check if index exists, else build it
if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${REFNAME}" ]; then
    echo "VARCAP:bwa_index for supplied reference exists, using it."
else
    echo "VARCAP:constructing bwa index."
    mkdir -p $PATH_BWA_075_INDEX/bwa_index_${REFNAME}
    cd $PATH_BWA_075_INDEX
    # bash $PATH_SCRIPTS/01_bwa_build_index.sh $REF_FA $PATH_BWA_075 $PATH_BWA_075_INDEX
    sbatch $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA} $PATH_BWA_075 $PATH_BWA_075_INDEX $PATH_LOGS >$PATH_LOGS/SLURM_B01.txt
    echo -e "qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA} $PATH_BWA_075 $PATH_BWA_075_INDEX" >>$PATH_LOGS/log.txt
fi

if [ -d "$REF_MAPPING" ]; then
  for EXNAME in $( ls "$REF_MAPPING" ); do
    EXNAME1=$( echo "$EXNAME" | sed 's/\.f.*a$//' )
    if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${EXNAME1}" ]; then
      echo "VARCAP: additional mapping bwa_index exists, using it."
    else
      # modify reference
      echo "Setup reference:$EXNAME"
      # ln -s "$REF_MAPPING/$EXNAME" "$REF_FA_MOD/${EXNAME}"
      # bash $PATH_PROJECTS_DATA/${PROJ_NAME}/A02_setup_reference.sh $REF_MAPPING/$EXNAME
      
      # construct bwa index
      echo "VARCAP:constructing additional mapping bwa index."
      mkdir -p $PATH_BWA_075_INDEX/bwa_index_${EXNAME1}
      cd $PATH_BWA_075_INDEX
      sbatch $PATH_SCRIPTS/A01_bwa_build_index.sh $REF_FA_MOD/${EXNAME} $PATH_BWA_075 ${PATH_BWA_075_INDEX} >>$PATH_LOGS/SLURM_B01.txt
    fi
  done
fi

echo -e "Bwa Index build." >>$PATH_PROJECTS_DATA/${PROJ_NAME}/logs/log.txt
