#!/bin/bash

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
if  ls "$PATH_LOGS"/SLURM_B02* 1> /dev/null 2>&1 ; then
  JOBLIST=$( cat $PATH_LOGS/SLURM_B02* | grep -o 'Submitted batch job [0-9]\+' | cut -d" " -f4 | tr "\n" ":" | sed 's/:$//' )
else
  echo "$( ls ${PATH_LOGS} )"
  echo -e "VARCAP_EXIT: No mapping job_id, in $PATH_LOGS/SLURM_B02* found.\nAssuming no mapping: Exiting.\n"
  exit 1
fi

# remove previous caller files
rm $PATH_LOFREQ_DATA/$OUTDIR/*
rm $PATH_LOFREQ2_DATA/$OUTDIR/*
rm $PATH_LOFREQ21_DATA/$OUTDIR/*
rm $PATH_VARSCAN_DATA/$OUTDIR/*
rm $PATH_FREEBAYES_DATA/$OUTDIR/*
rm -r $PATH_CORTEX_DATA/$OUTDIR/*

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
rm $PATH_LOGS/SLURM_C01_${i}.txt

### snp caller ###
# run samtools
echo "run samtools:"
mkdir -p $PATH_SAMTOOLS_DATA
mkdir -p $PATH_SAMTOOLS_DATA/$OUTDIR
cd $PATH_SAMTOOLS_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS
SAMEXE="$PATH_SCRIPTS/run_sam_bcftools.sh $REF_FNA $BAM_MOD $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $SAMEXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run gatk
echo "run gatk:"
mkdir -p $PATH_GATK_DATA
mkdir -p $PATH_GATK_DATA/$OUTDIR
cd $PATH_GATK_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS
GATEXE="$PATH_SCRIPTS/run_gatk.sh $REF_FA $BAM_MOD $PATH_GATK_DATA/$OUTDIR $PATH_GATK $PATH_PICARD $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $GATEXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run varscan
echo "run varscan:"
mkdir -p $PATH_VARSCAN_DATA
mkdir -p $PATH_VARSCAN_DATA/$OUTDIR
cd $PATH_VARSCAN_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS
VS2EXE="$PATH_SCRIPTS/run_varscan.sh $REF_FA $BAM_MOD $PATH_VARSCAN_DATA/$OUTDIR $PATH_VARSCAN $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $VS2EXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run lofreq
echo "run lofreq:"
mkdir -p $PATH_LOFREQ_DATA
mkdir -p $PATH_LOFREQ_DATA/$OUTDIR
cd $PATH_LOFREQ_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS
LOEXE="$PATH_SCRIPTS/run_lofreq.sh $REF_FA $BAM_MOD $PATH_LOFREQ_DATA/$OUTDIR $PATH_LOFREQ $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $LOEXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run lofreq2
echo "run lofreq2:"
mkdir -p $PATH_LOFREQ2_DATA
mkdir -p $PATH_LOFREQ2_DATA/$OUTDIR
cd $PATH_LOFREQ2_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS
LO2EXE="$PATH_SCRIPTS/run_lofreq2.sh $REF_FA $BAM_MOD $PATH_LOFREQ2_DATA/$OUTDIR $PATH_LOFREQ2 $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $LO2EXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run lofreq21
echo "run lofreq21:"
mkdir -p $PATH_LOFREQ21_DATA
mkdir -p $PATH_LOFREQ21_DATA/$OUTDIR
cd $PATH_LOFREQ21_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_lofreq21.sh $REF_FA $BAM_MOD $PATH_LOFREQ21_DATA/$OUTDIR $PATH_LOFREQ21 $PATH_SAMTOOLS
LO21EXE="$PATH_SCRIPTS/run_lofreq21.sh $REF_FA $BAM_MOD $PATH_LOFREQ21_DATA/$OUTDIR $PATH_LOFREQ21 $PATH_SAMTOOLS"
sbatch --dependency=afterok:${JOBLIST} $LO21EXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# run freebayes
echo "run freebayes:"
mkdir -p $PATH_FREEBAYES_DATA
mkdir -p $PATH_FREEBAYES_DATA/$OUTDIR
cd $PATH_FREEBAYES_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_freebayes.sh $REF_FA $BAM_MOD $PATH_FREEBAYES_DATA/$OUTDIR $PATH_FREEBAYES $PATH_SAMTOOLS
FBEXE="$PATH_SCRIPTS/run_freebayes.sh $REF_FA $BAM_MOD $PATH_FREEBAYES_DATA/$OUTDIR $PATH_FREEBAYES"
sbatch --dependency=afterok:${JOBLIST} $FBEXE >>$PATH_LOGS/SLURM_C01_${i}.txt

### sv caller ###
##run breakdancer
echo "run breakdancer:"
mkdir -p $PATH_BREAKD_DATA
mkdir -p $PATH_BREAKD_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH
BDREXE="$PATH_SCRIPTS/run_breakd_112_2013.sh $BAM_MOD $OUTDIR $PATH_BREAKD $PATH_BAM2CONF $PATH_SAMTOOLS $PATH_PICARD $PATH_BREAKD_DATA $PATH_GD_GRAPH_HISTOGRAM $READLENGTH"
sbatch --dependency=afterok:${JOBLIST} $BDREXE >>$PATH_LOGS/SLURM_C01_${i}.txt
##run pindel(024)
echo "run pindel:"
mkdir -p $PATH_PINDEL_DATA
mkdir -p $PATH_PINDEL_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE
PI1EXE="$PATH_SCRIPTS/run_pindel.sh $BAM_MOD $OUTDIR $PATH_PINDEL $PATH_SAM2PINDEL $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA $INSERT_SIZE"
sbatch --dependency=afterok:${JOBLIST} $PI1EXE >>$PATH_LOGS/SLURM_C01_${i}.txt
## run pindel_025
echo "run pindel_025:"
mkdir -p $PATH_PINDEL_DATA_025
mkdir -p $PATH_PINDEL_DATA_025/$OUTDIR
# sh $PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE
PI2EXE="$PATH_SCRIPTS/run_pindel_025.sh $BAM_MOD $OUTDIR $PATH_PINDEL_025 $PATH_SAM2PINDEL_025 $REF_FA $PATH_SAMTOOLS $PATH_PINDEL_DATA_025 $INSERT_SIZE"
sbatch --dependency=afterok:${JOBLIST} $PI2EXE >>$PATH_LOGS/SLURM_C01_${i}.txt
##run delly 0.7.2
echo "run delly 0.7.2:"
mkdir -p $PATH_DELLY_072_DATA
mkdir -p $PATH_DELLY_072_DATA/$OUTDIR
# sh $PATH_SCRIPTS/run_delly_072.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY_072 $PATH_DELLY_072_DATA $INSERT_SIZE $REF_FA_MOD
DLYEXE="$PATH_SCRIPTS/run_delly_072.sh $BAM_MOD $OUTDIR $REF_FA $PATH_DELLY_072 $PATH_DELLY_072_DATA $INSERT_SIZE $REF_FA_MOD"
sbatch --dependency=afterok:${JOBLIST} $DLYEXE >>$PATH_LOGS/SLURM_C01_${i}.txt


### comp assembly ###

# run cortex
echo "run cortex:"
# mkdir -p $PATH_CORTEX_DATA

##run cortex
cd $CURRENT_DIR
COREXE="$PATH_SCRIPTS/run_cortex_grid.sh $BAM_MOD $OUTDIR $i $BAM_BASE"
sbatch --dependency=afterok:${JOBLIST} $COREXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# bash $COREXE

# COREXE="$PATH_SCRIPTS/run_cortex_grid.sh  $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_1.fastq $PATH_CORTEX_DATA/$OUTDIR/${BAM_BASE}${i}_2.fastq $PATH_CORTEX_DATA/$OUTDIR/$BAM_BASE.$i $PATH_CORTEX_DATA \
# $REF_FA $PATH_VCFTOOLS $PATH_STAMPY $PATH_CORTEX_DIR $PATH_CORTEX $PATH_RUN_CALLS $PATH_SCRIPTS"
# qsub -l vf=5G -pe parallel 2 -sync n -hold_jid ${JOBLIST} -N waVCPCcor $COREXE >>$PATH_LOGS/SLURM_C01_${i}.txt
# rm ${BAM_BASE}${i}.rgroup.bam 

done
sleep 1


