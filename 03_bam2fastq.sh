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

# check for input
if [ "$#" -lt 1 ];then
  echo "Usage: bash bam2fastq.sh <file.bam>"
  echo "Please give absolute path to bam file."
  exit
fi

BAM_IN=$1
OUTDIR=$2
BAM_BASENAME=$( basename $BAM_IN .bam)

mkdir -p $OUTDIR
cd $OUTDIR

# create fastq files
java -jar -Xmx5g $PATH_PICARD/picard.jar SamToFastq INPUT=$BAM_IN FASTQ=${OUTDIR}/${BAM_BASENAME}_1.fastq SECOND_END_FASTQ=${OUTDIR}/${BAM_BASENAME}_2.fastq INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT

# check quality
# mkdir  -p fast_qc
# fastqc -o fast_qc/ ${BAM_BASENAME}_1.fastq
# fastqc -o fast_qc/ ${BAM_BASENAME}_2.fastq
# mv fast_qc $PATH_PROJECTS_DATA/$PROJ_NAME/raw/

echo "BAM2FASTQ finished" >>$PATH_PROJECTS_DATA/$PROJ_NAME/log.txt
