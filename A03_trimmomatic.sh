#!/bin/bash

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG


mkdir -p $PATH_PROJECTS_DATA/${PROJ_NAME}/filter/old_trimmomatic_logs
mv $PATH_PROJECTS_DATA/${PROJ_NAME}/filter/trimmomatic* $PATH_PROJECTS_DATA/${PROJ_NAME}/filter/old_trimmomatic_logs/ 2> /dev/null

# removes the illumina multiplex adapter from the reads and filters them according to the trimmomatic script
# cd ${PATH_PROJECTS_DATA}/${PROJ_NAME}/filter
echo "VARCAP: Trimmomatic started."
sbatch $PATH_SCRIPTS/run_trimmomatic.sh >$PATH_LOGS/SLURM_A03.txt

