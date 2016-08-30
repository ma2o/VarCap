#!/bin/bash

# SLURM
#SBATCH --job-name=VCbwaindex
#SBATCH --cpus-per-task=2
#SBATCH --mem=6000
#SBATCH --output=build_index-%A_%a.out
#SBATCH --error=build_index-%A_%a.err

REF=$1
REFNAME=$(basename $REF | sed 's/\.f.*$//')
PATH_BWA_2=$2
PATH_INDEX=$3
PATH_LOGS=$4

# write grid job_id to file for syncing
# echo "$JOB_ID started" >${PATH_LOGS}/B01_${JOB_ID}.txt

mkdir -p bwa_index_$REFNAME
# first step: index database
${PATH_BWA_2}/bwa index -p ${PATH_INDEX}/bwa_index_$REFNAME/$REFNAME $REF
# echo "$JOB_ID finished" >${PATH_LOGS}/B01_${JOB_ID}.txt

