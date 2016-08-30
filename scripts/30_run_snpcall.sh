#!/bin/bash
#$-q all.q@cube[ab]*

BAM=$1
OUTDIR=$2
REPEATS=$3

# snp caller
PATH_SAMTOOLS=$4
PATH_PICARD=$5
PATH_GATK=$6
PATH_VARSCAN=$7
# sv caller
PATH_BREAKD=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/breakdancer
PATH_PINDEL=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/pindel
PATH_DELLY=/scratch/zojer/projects/test_pipeline/mapping/recomb_calling/delly
# comp assembler
PATH_CORTEX=

REF_FNA=$8
REF_FA=$9
PATH_GATK_DATA=${10}
PATH_SAMTOOLS_DATA=${11}
PATH_VARSCAN_DATA=${12}

#run i times
for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "run: $i"
  #if repeat == 1, then do not modify name and run for once, else modify ending to 1...i and run ix
  if [ $REPEATS == 1 ]
    then
      BAM_MOD=$BAM
    else
      BAM_BASE=$(basename $BAM | sed 's/.\.bam$//')
      BAM_DIR=$(dirname $BAM)
      BAM_MOD=$BAM_DIR/${BAM_BASE}${i}.bam
  fi

echo "modified: "$BAM_MOD

### snp caller ###
#run samtools
echo "run samtools:"
mkdir $PATH_SAMTOOLS_DATA/$OUTDIR
cd $PATH_SAMTOOLS_DATA/$OUTDIR
sh $PATH_SAMTOOLS/run_snp_samtools.sh $REF_FNA $BAM_MOD $OUTDIR $PATH_SAMTOOLS
#run gatk
echo "run gatk:"
mkdir $PATH_GATK_DATA/$OUTDIR
cd $PATH_GATK_DATA/$OUTDIR
sh $PATH_GATK/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS
#run varscan
echo "run varscan:"
mkdir $PATH_VARSCAN_DATA/$OUTDIR
cd $PATH_VARSCAN_DATA/$OUTDIR
sh $PATH_VARSCAN/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS

### sv caller ###
#run breakdancer
mkdir $PATH_BREAKD/$OUTDIR
sh $PATH_BREAKD/run_breakd.sh $BAM_MOD $OUTDIR
#run pindel
mkdir $PATH_PINDEL/$OUTDIR
sh $PATH_PINDEL/run_pindel.sh $BAM_MOD $OUTDIR
#run delly
mkdir $PATH_DELLY/$OUTDIR
sh $PATH_DELLY/run_delly.sh $BAM_MOD $OUTDIR

### comp assembly ###
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

#check for empty files and delete them:
#check pindel dir for empty files and remove them
FILES_CHECK=$PATH_PINDEL/$OUTDIR/*
echo "Checking for empty files ..."
for f in $FILES_CHECK
do
  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
  if [ $FILESCK_SIZE == 0 ]
  then
    echo "remove "$f
    rm $f
  fi
done

#check delly dir for empty files and remove them
FILES_CHECK=$PATH_DELLY/$OUTDIR/*
echo "Checking for empty files ..."
for f in $FILES_CHECK
do
  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
  if [ $FILESCK_SIZE == 0 ]
  then
    echo "remove "$f
    rm $f
  fi
done

