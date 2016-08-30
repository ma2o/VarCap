#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# remove fq files from filter (except alt)
cd $PATH_PROJECTS_DATA/$PROJ_NAME/
find filter/ -name "*EMP_1*" -delete
find filter/ -name "*EMP_2*" -delete

# remove prinseq files
cd $PATH_PROJECTS_DATA/$PROJ_NAME/filter/
find prinseq_data/ -name "*.fastq" -delete

# remove subsample fastq files for mapping
cd $PATH_PROJECTS_DATA/$PROJ_NAME/mapper/
find subsample/ -name "*.fq*" -delete

# remove diverse files produced by calling
cd $PATH_PROJECTS_DATA/$PROJ_NAME/caller/
find . -name "*.ba?" -delete
find cortex/ -name "*.fastq" -delete
find cortex/ -name "*.ctx*" -delete


