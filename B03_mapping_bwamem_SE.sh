#!/bin/bash
#$-q all.q@cube[ab]*

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
if [ -f $PATH_LOGS/SGE_B01 ]; then
  JOBLIST=$( cat $PATH_LOGS/SGE_B01* $PATH_LOGS/SGE_A03* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
else
  JOBLIST=$( cat $PATH_LOGS/SGE_A03* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
fi

for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "bwa run: $i"
  OUTNAME_BWA=${BAM_NAME_BASE}_bwa
  echo "$OUTNAME_BWA"
  EXECOM="${PATH_SCRIPTS}/run_bwamem2bam_SE.sh $PATH_PROJECTS_DATA/$PROJ_NAME $i $SYSTEM_CONFIG"
  qsub -l vf=4G -pe parallel 2 ${EXECOM}
  # qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPMbwa ${EXECOM} >$PATH_LOGS/SGE_B02_${i}.txt
  # qsub -l vf=4G -pe parallel 2 ${PATH_SCRIPTS}/run_bwamem2bam.sh ${REF_FA}_map.fasta ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_1.fq ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_2.fq ${OUTNAME_BWA}_v${i} $PATH_BWA_075 $PATH_BWA_DATA $PATH_BWA_075_INDEX $PATH_PICARD $PATH_SAMTOOLS $PATH_LOGS >$PATH_LOGS/SGE_B02_${i}.txt
done


