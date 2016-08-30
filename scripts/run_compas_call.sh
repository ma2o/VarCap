#!/bin/bash
#$-q all.q@cube[ab]*
BAM=$1
OUTDIR=$2
REPEATS=$3

PATH_PICARD=/home/apps/picard-tools-1.62
PATH_SAMTOOLS=/home/apps/samtools-0.1.18
PATH_CORTEX=/scratch/zojer/projects/test_pipeline/mapping/comp_assembly/cortex

#run i times
for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "run: $i"
  #if repeat == 1, then do not modify name and run for once, else modify ending to 1...i and run ix
  if [ $REPEATS == 1 ]
    then
      BAM_MOD=$BAM
      BAM_BASE=$(basename $BAM | sed 's/\..am$//')
      i=1
    else
      BAM_BASE=$(basename $BAM | sed 's/.\..am$//')
      BAM_DIR=$(dirname $BAM)
      BAM_MOD=$BAM_DIR/${BAM_BASE}${i}.bam
  fi

echo "modified: "$BAM_MOD

#run cortex pipeline
mkdir $PATH_CORTEX/$OUTDIR/$BAM_BASE.$i
cd $PATH_CORTEX/$OUTDIR

# replace/add read group header
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/AddOrReplaceReadGroups.jar I=$BAM_MOD O=${BAM_BASE}${i}.rgroup.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false
#convert bam file to fastq reads
java -jar -Xmx3g $PATH_PICARD/SamToFastq.jar INPUT=${BAM_BASE}${i}.rgroup.bam FASTQ=${BAM_BASE}${i}_1.fastq SECOND_END_FASTQ=${BAM_BASE}${i}_2.fastq INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
#run cortex
sh $PATH_CORTEX/run_cortex.sh  $PATH_CORTEX/$OUTDIR/${BAM_BASE}${i}_1.fastq $PATH_CORTEX/$OUTDIR/${BAM_BASE}${i}_2.fastq $PATH_CORTEX/$OUTDIR/$BAM_BASE.$i
rm ${BAM_BASE}${i}.rgroup.bam ${BAM_BASE}${i}_1.fastq ${BAM_BASE}${i}_2.fastq

done
