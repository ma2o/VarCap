#!/bin/bash

# SLURM
#SBATCH --job-name=VCfastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --output=raw/fastqc-%A_%a.out
#SBATCH --error=raw/fastqc-%A_%a.err

. /etc/profile

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

INFILE_RAW=$1
if [[ -z $INFILE_RAW ]]; then
  # if no file is supplied, run quality checks for raw data folder
  INFILE_RAW=($( find  $PATH_PROJECTS_DATA/${PROJ_NAME}/raw/* | grep -E 'fq$|fastq$|gz$|bam$'))
  for line in "${INFILE_RAW[@]}"; do
    # run fastqc
    OUTDIR_QC=fastqc_raw
    cd $PATH_PROJECTS_DATA/${PROJ_NAME}/raw
    mkdir -p $OUTDIR_QC
    $PATH_FASTQC/fastqc -o $OUTDIR_QC/ $INFILE_RAW
    # print log messages
    echo -e "FastQC raw finished: $INFILE_RAW" >>$PATH_LOGS/log.txt
    echo -e "$PATH_FASTQC/fastqc -o $OUTDIR_QC/ $INFILE_RAW" >>$PATH_LOGS/log.txt
  done
else
  # run fastqc
  OUTDIR_QC=fastqc_raw
  cd $PATH_PROJECTS_DATA/${PROJ_NAME}/raw
  mkdir -p $OUTDIR_QC
  $PATH_FASTQC/fastqc -o $OUTDIR_QC/ $INFILE_RAW

  # print log messages
  echo -e "FastQC raw finished: $INFILE_RAW" >>$PATH_LOGS/log.txt
  echo -e "$PATH_FASTQC/fastqc -o $OUTDIR_QC/ $INFILE_RAW" >>$PATH_LOGS/log.txt
fi

