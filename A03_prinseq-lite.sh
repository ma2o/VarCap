#!/bin/bash

# Script filters paired fastq reads.

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

### 11 check read names, and modify if they are in invalid format
# NEW_SUFFIX_1="/1"
# NEW_SUFFIX_2="/2"
OUTSUFF1=alt_1
OUTSUFF2=alt_2

cd $PATH_PROJECTS_DATA/${PROJ_NAME}/filter
# echo "Correct readnames."
# perl $PATH_SCRIPTS/correct_readname.pl ${PROJ_NAME}_ad_1.fastq $NEW_SUFFIX_1 >${PROJ_NAME}_1.fastq
# perl $PATH_SCRIPTS/correct_readname.pl ${PROJ_NAME}_ad_2.fastq $NEW_SUFFIX_2 >${PROJ_NAME}_2.fastq

### 12 quality filter using prinseq-lite
echo "Quality filtering."
FILENR=$( ls $PATH_PROJECTS_DATA/$PROJ_NAME/raw | grep -E 'fastq$|fq$' | wc -l )
COUNT=1
for file in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/raw | grep -E 'fastq$|fq$' ); do
  echo "$COUNT of $FILENR:$file"
  sbatch $PATH_SCRIPTS/13_qfilter_comp.sh $PATH_PROJECTS_DATA/$PROJ_NAME/raw/$file $READS1_TRIM $READS1_TRIMF $QUAL_WINDOW_1 $MIN_LENGTH_1 $PATH_PRINSEQ $PATH_FASTQC $OUTSUFF1
  if [[ $COUNT -eq 2 ]]; then
    echo "2 of $FILENR:$file"
    sbatch $PATH_SCRIPTS/13_qfilter_comp.sh $PATH_PROJECTS_DATA/$PROJ_NAME/raw/$file $READS2_TRIM $READS2_TRIMF $QUAL_WINDOW_2 $MIN_LENGTH_2 $PATH_PRINSEQ $PATH_FASTQC $OUTSUFF2
  fi
  let COUNT=$COUNT+1
done

### 13 to pairs
# echo "Extract read pairs."
# sh $PATH_SCRIPTS/filter/read2pairs/red2readpairs.sh ${PROJ_NAME}_1_filter.fastq ${PROJ_NAME}_2_filter.fastq

### 14 rename/remove files, first check if read pairs file exist, then proceed with file modification

# if [ -f *_2_filter_rp.fastq ];
# then
  # rm ${PROJ_NAME}_1.fastq ${PROJ_NAME}_2.fastq *_filter.fastq
  # mv *_1_filter_rp.fastq ${PROJ_NAME}_alt_1.fastq
  # mv *_2_filter_rp.fastq ${PROJ_NAME}_alt_2.fastq
  # rm prinseq_data/*.fastq
# fi

echo "Filter step finished. - You may proceed with 20 subsampling/mapping."
