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

# check if we have enough reads for subsampling, else take everything we have
BS=$( cat $PATH_PROJECTS_DATA/${PROJ_NAME}/info.txt | grep -e 'BothS' | cut -d' ' -f2 )
if [ "$SUBSAMPLE_SIZE_ALT" -gt "$BS" ]; then
    sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$BS"'/' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available."
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available." >>$PATH_PROJECTS_DATA/${PROJ_NAME}/log.txt
    SUBSAMPLE_SIZE_ALT=$BS
fi

REFNAME=$(basename $REF_FA | sed 's/\.f.*$//')

#subsample simulated reads (path to reference/variant reads set in subscript)
if [ ! -d $PATH_SUBSAMPLE_DATA ];
  then
  mkdir -p $PATH_SUBSAMPLE_DATA
fi
cd $PATH_SUBSAMPLE_DATA

bash ${PATH_SCRIPTS}/22_sample_subsets.sh $SUBSAMPLE_SIZE_ALT $BAM_NAME_BASE $REPEATS $PATH_ALT_READS1 $PATH_ALT_READS2 $PATH_SCRIPTS

### bwa ---
## prepare mapper (index)
# create bwa data folder
mkdir -p $PATH_BWA_DATA

  # check if index exists, else build it
  if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${REFNAME}" ]; then
    echo "VARCAP:bwa_index exists, using it."
  else
    echo "VARCAP:constructing bwa index."
    mkdir -p $PATH_BWA_075_INDEX/bwa_index_${REFNAME}
    cd $PATH_BWA_075_INDEX
    bash $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA}_map.fasta $PATH_BWA_075 $PATH_BWA_075_INDEX
    # qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA}_map.fasta $PATH_BWA_075 $PATH_BWA_075_INDEX
  fi

# run bwa
cd $PATH_BWA_DATA

# get job_id list of submitted SGI jobs for syncing
JOBLIST=$( cat $PATH_LOGS/SLURM_B01* $PATH_LOGS/SLURM_A03* | grep -o 'Your job [0-9]\+' | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )

for ((  i = 1 ; i <= $REPEATS;  i++  ))
do
  echo "bwa run: $i"
  OUTNAME_BWA=${BAM_NAME_BASE}_bwa
  echo "$OUTNAME_BWA"
  $EXECOM="${PATH_SCRIPTS}/run_bwamem2bam.sh ${REF_FA}_map.fasta ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_1.fq ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_2.fq ${OUTNAME_BWA}_v${i} $PATH_BWA_075 $PATH_BWA_DATA $PATH_BWA_075_INDEX $PATH_PICARD $PATH_SAMTOOLS $PATH_LOGS >$PATH_LOGS/SGE_B02_${i}.txt"
  qsub -l vf=4G -pe parallel 2 -sync y -hold_jid ${JOBLIST} -N waitVCP 'sleep 2' ${EXECOM}
  # qsub -l vf=4G -pe parallel 2 ${PATH_SCRIPTS}/run_bwamem2bam.sh ${REF_FA}_map.fasta ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_1.fq ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${i}_2.fq ${OUTNAME_BWA}_v${i} $PATH_BWA_075 $PATH_BWA_DATA $PATH_BWA_075_INDEX $PATH_PICARD $PATH_SAMTOOLS $PATH_LOGS >$PATH_LOGS/SGE_B02_${i}.txt
done


