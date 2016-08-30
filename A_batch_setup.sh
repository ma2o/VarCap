#!/bin/bash
#$-q all.q@cube[ab]*

INDIR=$1
REGEX_FILTER=$2
VARCAP_DIR=/scratch/zojer/projects/varcap_2.1
PATH_PROJECTS_DATA=/scratch/zojer/projects/varcap_data2
REF=/proj/evochlamy/ref_upd/NC_005861_PAC_03_mito.fasta

for file in $( ls $INDIR | grep -E '^[0-9]' | grep -E "${REGEX_FILTER}" ); do 
  echo $file;
  NAME=$( echo $file | sed 's/_[0-9,A-Z]*_[0-9,A-Z]*.bam$//');
  echo $NAME;
  # create projects
  cd $VARCAP_DIR
  bash 00_setup_project.sh $NAME $PATH_PROJECTS_DATA $REF
  cd $PATH_PROJECTS_DATA/$NAME
  bash 01_setup_reference.sh
  bash 02_setup_index.sh
  bash 03_bam2fastq.sh ${INDIR}/$file
  # bash 1101_adapter_trimmomatic.sh
  
done


