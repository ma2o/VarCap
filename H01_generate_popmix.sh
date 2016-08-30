#!/bin/bash
#$-q all.q@cube[ab]*
#$ -o log.fastqc.raw
#$ -e log.fastqc.raw

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

OUTDIR=$1
if [[ -z $OUTDIR ]]; then
  OUTDIR=mix
fi

mkdir -p $PATH_PROJECTS_DATA/${PROJ_NAME}/mix
cd $PATH_PROJECTS_DATA/${PROJ_NAME}/mix

bash ${PATH_SCRIPTS}/23_sample_subsets.sh $SUBSAMPLE_SIZE_B $SUBSAMPLE_SIZE_A $BAM_NAME_BASE $OUTDIR $PATH_A_READS1 $PATH_A_READS2 $PATH_B_READS1 $PATH_B_READS2 $PATH_SCRIPTS
