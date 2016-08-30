#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

if [ "$#" -lt 5 ];then
  echo "Converts vcf file to ccf (column only) file format"
  echo "Usage: bash 407_vcf2ccf.sh <project_id> <experiment_id> <sequencing_id> <experiment> <replicate> <selection> <selection_value> <timepoint>"
  exit
fi

PROJ_ID=$1
EXP_ID=$2
SEQ_ID=$3
EXPERIMENT=$4
REPLICATE=$5
SELECTION=$6
SELECTION_VALUE=$7
TIMEPOINT=$8

# FILE_NAME_BASE=$( basename $FILE_NAME | sed 's/\.vcf$//' )

cd $PATH_PROJECTS_DATA/$PROJ_NAME

### filter vcfs: apply cutoff parameters in form of absolute read count/ relative percentage coverage e.g. 16 8 (lower limit)
echo "VARCAP: convert vcf to ccf"
# cd $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw
# bash $PATH_SCRIPTS/reformat_vcf2ccf.sh $EXP_NAME $FILE_NAME $PATH_SCRIPTS >${FILE_NAME_BASE}.ccf

cd $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs
COUNTER=1
# convert all vcf files in folder annotator/snpeff to ccf
for file in $PATH_PROJECTS_DATA/$PROJ_NAME/annotator/snpeff/*.vcf
do
  echo "$file"
  FILE_BASE=$( basename $file)
  FILE_NAME_BASE=${FILE_BASE%.vcf}
  bash $PATH_SCRIPTS/reformat_vcf2ccf.sh $PROJ_ID $EXP_ID $SEQ_ID $EXPERIMENT $REPLICATE $SELECTION $SELECTION_VALUE $TIMEPOINT $FILE_NAME_BASE $COUNTER $file $PATH_SCRIPTS >${FILE_NAME_BASE}.ccf
  let COUNTER=COUNTER+1
done

