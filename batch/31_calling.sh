#!/bin/bash

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

RC_RATIO=$(((100 * ${SUBSAMPLE_SIZE_ALT}) / (${SUBSAMPLE_SIZE_REF} + ${SUBSAMPLE_SIZE_ALT}) | bc))
BAM=$PATH_BWA_DATA/${BAM_NAME_BASE}_bwa_${RC_RATIO}_v1.bam
# BAM=$PATH_BAM_BWA_FOR_CALLER
OUTDIR=$OUTPUT_DIR

# run i times
for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "run: $i"
  #if repeat == 1, then do not modify name and run for once, else modify ending to 1...i and run ix
  if [ $REPEATS == 1 ]
    then
      BAM_MOD=$BAM
      BAM_BASE=$(basename $BAM | sed 's/.\.bam$//')
    else
      BAM_BASE=$(basename $BAM | sed 's/.\.bam$//')
      BAM_DIR=$(dirname $BAM)
      BAM_MOD=$BAM_DIR/${BAM_BASE}${i}.bam
  fi

echo "modified: "$BAM_MOD

### snp caller ###
#run samtools
# echo "run samtools:"
mkdir -p $PATH_SAMTOOLS_DATA
mkdir -p $PATH_SAMTOOLS_DATA/$OUTDIR
cd $PATH_SAMTOOLS_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS
# qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS
#run gatk
# echo "run gatk:"
mkdir -p $PATH_GATK_DATA
mkdir -p $PATH_GATK_DATA/$OUTDIR
cd $PATH_GATK_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS
# qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS
#run varscan
echo "run varscan:"
if [ ! -d $PATH_VARSCAN_DATA ];then
  mkdir -p $PATH_VARSCAN_DATA
fi
mkdir -p $PATH_VARSCAN_DATA/$OUTDIR
cd $PATH_VARSCAN_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS
qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS
#run lofreq
echo "run lofreq:"
if [ ! -d $PATH_LOFREQ_DATA ];then
  mkdir -p $PATH_LOFREQ_DATA
fi
mkdir -p $PATH_LOFREQ_DATA/$OUTDIR
cd $PATH_LOFREQ_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS
#run lofreq2
echo "run lofreq2:"
if [ ! -d $PATH_LOFREQ2_DATA ];then
  mkdir -p $PATH_LOFREQ2_DATA
fi
mkdir -p $PATH_LOFREQ2_DATA/$OUTDIR
cd $PATH_LOFREQ2_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS


### sv caller ###
##run breakdancer
echo "run breakdancer:"
if [ ! -d $PATH_BREAKD_DATA ];then
  mkdir -p $PATH_BREAKD_DATA
fi
mkdir -p $PATH_BREAKD_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH
##run pindel(024)
echo "run pindel:"
if [ ! -d $PATH_PINDEL_DATA ];then
  mkdir -p $PATH_PINDEL_DATA
fi
mkdir -p $PATH_PINDEL_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE
## run pindel_025
echo "run pindel_025:"
if [ ! -d $PATH_PINDEL_DATA_025 ];then
  mkdir -p $PATH_PINDEL_DATA_025
fi
mkdir -p $PATH_PINDEL_DATA_025/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE
##run delly
echo "run delly:"
if [ ! -d $PATH_DELLY_DATA ];then
  mkdir -p $PATH_DELLY_DATA
fi
mkdir -p $PATH_DELLY_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_delly.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY $PATH_DELLY_DATA $INSERT_SIZE $REF_FA_MOD
qsub -l vf=4G -pe parallel 2 $PATH_SCRIPTS/run_delly.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY $PATH_DELLY_DATA $INSERT_SIZE $REF_FA_MOD

### comp assembly ###
##run cortex pipeline
echo "run cortex:"
if [ ! -d $PATH_CORTEX_DATA ];then
  mkdir -p $PATH_CORTEX_DATA
fi
mkdir -p $PATH_CORTEX_DATA/$OUTDIR
mkdir -p $PATH_CORTEX_DATA/$OUTDIR/$BAM_BASE.$i
cd $PATH_CORTEX_DATA/$OUTDIR
## replace/add read group header
# java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/AddOrReplaceReadGroups.jar I=$BAM_MOD O=${BAM_BASE}${i}.rgroup.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false
##convert bam file to fastq reads
# java -jar -Xmx3g $PATH_PICARD/SamToFastq.jar INPUT=${BAM_BASE}${i}.rgroup.bam FASTQ=${BAM_BASE}${i}_1.fastq SECOND_END_FASTQ=${BAM_BASE}${i}_2.fastq INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
##run cortex
# qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/run_cortex.sh  $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_1.fastq $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_2.fastq $PATH_CORTEX_DATA/$OUTDIR/$BAM_BASE.$i $PATH_CORTEX_DATA \
# $REF_FA $PATH_VCFTOOLS $PATH_STAMPY $PATH_CORTEX_DIR $PATH_CORTEX $PATH_RUN_CALLS
# rm ${BAM_BASE}${i}.rgroup.bam 

done

##check for empty files and delete them:
##check pindel dir for empty files and remove them
#FILES_CHECK=$PATH_PINDEL/$OUTDIR/*
#echo "Checking for empty files ..."
#for f in $FILES_CHECK
#do
#  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
#  if [ $FILESCK_SIZE == 0 ]
#  then
#    echo "remove "$f
#    rm $f
#  fi
#done

###check delly dir for empty files and remove them
#FILES_CHECK=$PATH_DELLY/$OUTDIR/*
#echo "Checking for empty files ..."
#for f in $FILES_CHECK
#do
#  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
#  if [ $FILESCK_SIZE == 0 ]
#  then
#    echo "remove "$f
#    rm $f
#  fi
#done

