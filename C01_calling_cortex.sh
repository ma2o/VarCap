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

BAM=$PATH_BWA_DATA/${BAM_NAME_BASE}_bwa_1.bam
# BAM=$PATH_BAM_BWA_FOR_CALLER
OUTDIR=$OUTPUT_DIR

# check for mapper output/bam file and create joblist of job_id for syncing of SGE
if  ls "$PATH_LOGS"/SGE_B02* 1> /dev/null 2>&1 ; then
  JOBLIST=$( cat $PATH_LOGS/SGE_B02* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
else
  echo "$( ls ${PATH_LOGS} )"
  echo -e "VARCAP_EXIT: No mapping job_id, in $PATH_LOGS/SGE_B02* found.\nAssuming no mapping: Exiting.\n"
  exit 1
fi


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

# remove old job logs
rm $PATH_LOGS/SGE_C01_${i}.txt

### snp caller ###
#run samtools
# echo "run samtools:"
mkdir -p $PATH_SAMTOOLS_DATA
mkdir -p $PATH_SAMTOOLS_DATA/$OUTDIR
cd $PATH_SAMTOOLS_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS
SAMEXE="$PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS"
# qsub -l vf=5G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCsam $SAMEXE >>$PATH_LOGS/SGE_C01_${i}.txt
#run gatk
# echo "run gatk:"
mkdir -p $PATH_GATK_DATA
mkdir -p $PATH_GATK_DATA/$OUTDIR
cd $PATH_GATK_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS
GATEXE="$PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS"
# qsub -l vf=5G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCgat $GATEXE >>$PATH_LOGS/SGE_C01_${i}.txt
#run varscan
# echo "run varscan:"
if [ ! -d $PATH_VARSCAN_DATA ];then
  mkdir -p $PATH_VARSCAN_DATA
fi
mkdir -p $PATH_VARSCAN_DATA/$OUTDIR
cd $PATH_VARSCAN_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS
VS2EXE="$PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS"
# qsub -l vf=5G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCvs2 $VS2EXE >>$PATH_LOGS/SGE_C01_${i}.txt
#run lofreq
# echo "run lofreq:"
if [ ! -d $PATH_LOFREQ_DATA ];then
  mkdir -p $PATH_LOFREQ_DATA
fi
mkdir -p $PATH_LOFREQ_DATA/$OUTDIR
cd $PATH_LOFREQ_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS
LO1EXE="$PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPClo1 $LO1EXE >>$PATH_LOGS/SGE_C01_${i}.txt
#run lofreq2
# echo "run lofreq2:"
if [ ! -d $PATH_LOFREQ2_DATA ];then
  mkdir -p $PATH_LOFREQ2_DATA
fi
mkdir -p $PATH_LOFREQ2_DATA/$OUTDIR
cd $PATH_LOFREQ2_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS
LO2EXE="$PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPClo2 $LO2EXE >>$PATH_LOGS/SGE_C01_${i}.txt


### sv caller ###
##run breakdancer
# echo "run breakdancer:"
if [ ! -d $PATH_BREAKD_DATA ];then
  mkdir -p $PATH_BREAKD_DATA
fi
mkdir -p $PATH_BREAKD_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH
BDREXE="$PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCbdr $BDREXE >>$PATH_LOGS/SGE_C01_${i}.txt
##run pindel(024)
# echo "run pindel:"
if [ ! -d $PATH_PINDEL_DATA ];then
  mkdir -p $PATH_PINDEL_DATA
fi
mkdir -p $PATH_PINDEL_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE
PI1EXE="$PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCpi1 $PI1EXE >>$PATH_LOGS/SGE_C01_${i}.txt
## run pindel_025
# echo "run pindel_025:"
if [ ! -d $PATH_PINDEL_DATA_025 ];then
  mkdir -p $PATH_PINDEL_DATA_025
fi
mkdir -p $PATH_PINDEL_DATA_025/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE
PI2EXE="$PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCpi2 $PI2EXE >>$PATH_LOGS/SGE_C01_${i}.txt
##run delly 0.0.9
# echo "run delly 0.0.9:"
if [ ! -d $PATH_DELLY_DATA ];then
  mkdir -p $PATH_DELLY_DATA
fi
mkdir -p $PATH_DELLY_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_delly.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY $PATH_DELLY_DATA $INSERT_SIZE $REF_FA_MOD
DLYEXE="$PATH_SCRIPTS/run_delly.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY $PATH_DELLY_DATA $INSERT_SIZE $REF_FA_MOD"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCdly $DLYEXE >>$PATH_LOGS/SGE_C01_${i}.txt
##run delly 0.7.2
# echo "run delly 0.7.2:"
mkdir -p $PATH_DELLY_072_DATA
mkdir -p $PATH_DELLY_072_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_delly_072.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY_072 $PATH_DELLY_072_DATA $INSERT_SIZE $REF_FA_MOD
DLYEXE="$PATH_SCRIPTS/run_delly_072.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY_072 $PATH_DELLY_072_DATA $INSERT_SIZE $REF_FA_MOD"
# qsub -l vf=4G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCdly $DLYEXE >>$PATH_LOGS/SGE_C01_${i}.txt


### comp assembly ###

##run cortex pipeline

# prepare cortex index
# mkdir -p $PATH_CORTEX_DATA
# cd $PATH_CORTEX_DATA
# rm $PATH_CORTEX_DATA/*
# rm $PATH_CORTEX_DATA/ref/*
# bash $PATH_SCRIPTS/00_cortex_prepare.sh $PATH_CORTEX_DATA $REF_FA $PATH_VCFTOOLS $PATH_STAMPY $PATH_CORTEX_DIR $PATH_CORTEX $PATH_RUN_CALLS

# run cortex
echo "run cortex:"
# mkdir -p $PATH_CORTEX_DATA

# mkdir -p $PATH_CORTEX_DATA/$OUTDIR
# mkdir -p $PATH_CORTEX_DATA/$OUTDIR/$BAM_BASE.$i
# cd $PATH_CORTEX_DATA/$OUTDIR
## replace/add read group header
# java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/AddOrReplaceReadGroups.jar I=$BAM_MOD O=${BAM_BASE}${i}.rgroup.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false
##convert bam file to fastq reads
# java -jar -Xmx3g $PATH_PICARD/SamToFastq.jar INPUT=${BAM_BASE}${i}.rgroup.bam FASTQ=${BAM_BASE}${i}_1.fastq SECOND_END_FASTQ=${BAM_BASE}${i}_2.fastq INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
##run cortex
cd $CURRENT_DIR
COREXE="$PATH_SCRIPTS/run_cortex_grid.sh $BAM_MOD $OUTDIR $i $BAM_BASE"
qsub -l vf=6G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCcor $COREXE >>$PATH_LOGS/SGE_C01_${i}.txt
# bash $COREXE

# COREXE="$PATH_SCRIPTS/run_cortex_grid.sh  $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_1.fastq $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_2.fastq $PATH_CORTEX_DATA/$OUTDIR/$BAM_BASE.$i $PATH_CORTEX_DATA \
# $REF_FA $PATH_VCFTOOLS $PATH_STAMPY $PATH_CORTEX_DIR $PATH_CORTEX $PATH_RUN_CALLS $PATH_SCRIPTS"
# qsub -l vf=5G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCcor $COREXE >>$PATH_LOGS/SGE_C01_${i}.txt
# rm ${BAM_BASE}${i}.rgroup.bam 

done
sleep 1
## check for empty files and delete them:
## check pindel dir for empty files and remove them
# FILES_CHECK=$PATH_PINDEL/$OUTDIR/*
# echo "Checking for empty files ..."
# for f in $FILES_CHECK
# do
#  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
#  if [ $FILESCK_SIZE == 0 ]
#  then
#    echo "remove "$f
#    rm $f
#  fi
# done

### check delly dir for empty files and remove them
# FILES_CHECK=$PATH_DELLY/$OUTDIR/*
# echo "Checking for empty files ..."
# for f in $FILES_CHECK
# do
#  FILESCK_SIZE=$(du $f | sed 's/^\(.\).*$/\1/')
#  if [ $FILESCK_SIZE == 0 ]
#  then
#    echo "remove "$f
#    rm $f
#  fi
# done

