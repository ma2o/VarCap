#!/bin/bash

REGEX_FILTER=$1

SAMDIR=/home/apps/samtools-0.1.18
ACT_DIR=$( pwd )
for file in $( ls -d */ | grep -E "${REGEX_FILTER}" ); do
  echo $file;
  cd $file/mapper/bwa
  BAM_NAME=$( ls | grep -e '_all_.*bam' | sed 's/\.bam//' )
  mv ${BAM_NAME}.bam ${BAM_NAME}.map.bam
  mv ${BAM_NAME}.bai ${BAM_NAME}.map.bai
  $SAMDIR/samtools view -b ${BAM_NAME}.map.bam NC_000117.1 >${BAM_NAME}.bam
  $SAMDIR/samtools index ${BAM_NAME}.bam ${BAM_NAME}.bai
  cd $ACT_DIR
done

