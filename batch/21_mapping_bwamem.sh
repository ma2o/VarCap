#!/bin/bash

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

REFNAME=$(basename $REF_FA | sed 's/\.f.*$//')
RC_RATIO=$(((100 * ${SUBSAMPLE_SIZE_ALT}) / (${SUBSAMPLE_SIZE_REF} + ${SUBSAMPLE_SIZE_ALT}) | bc))

#subsample simulated reads (path to reference/variant reads set in subscript)
if [ ! -d $PATH_SUBSAMPLE_DATA ];
  then
  mkdir -p $PATH_SUBSAMPLE_DATA
fi
cd $PATH_SUBSAMPLE_DATA

bash ${PATH_SCRIPTS}/22_sample_subsets.sh $SUBSAMPLE_SIZE_REF $SUBSAMPLE_SIZE_ALT $BAM_NAME_BASE $REPEATS $PATH_ALT_READS1 $PATH_ALT_READS2 $PATH_REF_READS1 $PATH_REF_READS2 $PATH_SCRIPTS

### bwa ---
## prepare mapper (index)
# create bwa data folder
mkdir -p $PATH_BWA_DATA

  # check if index exists, else build it
  if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${REFNAME}" ]; then
    echo "VARCAP:bwa_index exists, using it."
  else
    echo "VARCAP:constructing bwa index."
    mkdir -p $PATH_BWA_075_INDEX/bwa_index_$REF_FA
    cd $PATH_BWA_075_INDEX
    # bash $PATH_SCRIPTS/01_bwa_build_index.sh $REF_FA $PATH_BWA_075 $PATH_BWA_075_INDEX
    qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/A01_bwa_build_index.sh $REF_FA $PATH_BWA_075 $PATH_BWA_075_INDEX
  fi

# run bwa
cd $PATH_BWA_DATA

for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "bwa run: $i"
  OUTNAME_BWA=${BAM_NAME_BASE}_bwa
  echo $OUTNAME_BWA
  qsub -l vf=4G -pe parallel 2 ${PATH_SCRIPTS}/run_bwamem2bam.sh $REF_FA ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_${RC_RATIO}_r${i}_1.fq ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_${RC_RATIO}_r${i}_2.fq ${OUTNAME_BWA}_${RC_RATIO}_v${i} $PATH_BWA_075 $PATH_BWA_DATA $PATH_BWA_075_INDEX $PATH_PICARD $PATH_SAMTOOLS
done


