#!/bin/bash

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG


### bwa ---
# 1. subsample
# 2. check else prepare index
# 3. run bwamem

# create bwa data folder
mkdir -p $PATH_BWA_DATA
cd $PATH_BWA_DATA

# get job_id list of submitted SGI jobs for syncing
# if [ -f $PATH_LOGS/SGE_B01 ]; then
#   JOBLIST=$( cat $PATH_LOGS/SGE_B01* $PATH_LOGS/SGE_A03* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
# else
#   JOBLIST=$( cat $PATH_LOGS/SGE_A03* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
# fi

# get job_id list of submitted SLURM jobs for syncing
# SLURM output: Submitted batch job 98765
if [ -f $PATH_LOGS/SLURM_B01 ]; then
  JOBLIST=$( cat $PATH_LOGS/SLURM_B01* $PATH_LOGS/SLURM_A03* | grep -o 'Submitted batch job [0-9]\+' | cut -d" " -f4 | tr "\n" ":" | sed 's/:$//' )
else
  JOBLIST=$( cat $PATH_LOGS/SLURM_A03* | grep -o 'Submitted batch job [0-9]\+' | cut -d" " -f4 | tr "\n" ":" | sed 's/:$//' )
fi

for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "bwa run: $i"
  OUTNAME_BWA=${BAM_NAME_BASE}_bwa
  echo "$OUTNAME_BWA"
  EXECOM="${PATH_SCRIPTS}/run_bwamem2bam2.sh $PATH_PROJECTS_DATA/$PROJ_NAME $i $SYSTEM_CONFIG"
  # bash ${EXECOM}
  # --dependency=afterok:<jobID_A:jobID_C:jobID_D>
  echo ${JOBLIST}
  sbatch --dependency=afterok:${JOBLIST} ${EXECOM} >$PATH_LOGS/SLURM_B02_${i}.txt
done


