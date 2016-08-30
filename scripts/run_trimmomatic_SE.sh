#!/bin/bash
#$-q all.q@cube[ab]*
#$ -o filter
#$ -e filter

# read config file(variant.config) within the same directory
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# check if files are in bam format and convert to fastq if necessary
INFILE_BAM=$( find  $PATH_PROJECTS_DATA/$PROJ_NAME/raw/* | grep -E 'fq$|fastq$|gz$|bam$' | head -n1 )
OUTDIR_TEMP=$TMPDIR/varcap/${PROJ_NAME}
if [[ $INFILE_BAM == *bam ]]; then
  mkdir -p $OUTDIR_TEMP
  bash $PATH_VARCAP/03_bam2fastq.sh $INFILE_BAM $OUTDIR_TEMP
  echo -e "$PATH_VARCAP/03_bam2fastq.sh $INFILE_BAM $OUTDIR_TEMP"
  # create links
  BAM_BASENAME=$( basename $BAM_IN .bam)
  for file in $( find $OUTDIR_TEMP/* | grep -E 'fastq$|gz$|' | sed 's/.*\///' ); do
      cd $PATH_PROJECTS_DATA/$PROJ_NAME/raw
      ln -s $OUTDIR_TEMP/$file $file
  done
fi

### trimmomatic and fastqc
OUTDIR_QC=fastqc
TAILCROP=$(( $READLENGTH - $READS1_TRIM ))
# removes the illumina multiplex adapter from the reads, per default the multiplex identifier is set to NNNNNN
# however you can provide the exact sequence within the command line


mkdir -p $PATH_PROJECTS_DATA/${PROJ_NAME}/filter
cd $PATH_PROJECTS_DATA/${PROJ_NAME}/filter

READS1=$PATH_PROJECTS_DATA/${PROJ_NAME}/raw/$( ls  $PATH_PROJECTS_DATA/${PROJ_NAME}/raw | grep -e '1\.f.*q.*$' )
READS2=$PATH_PROJECTS_DATA/${PROJ_NAME}/raw/$( ls  $PATH_PROJECTS_DATA/${PROJ_NAME}/raw | grep -e '2\.f.*q.*$' )

echo -e "READS1:$READS1" >>$PATH_LOGS/log.txt
echo -e "READS2:$READS2" >>$PATH_LOGS/log.txt


if [[ $READS1_TRIM -eq 0 ]]; then
  echo "Run Trimmomatic no readtrim."
  java -jar $PATH_TRIMMOMATIC/trimmomatic-0.32.jar PE -threads 2 -phred33  ${READS1} ${READS2} ${PROJ_NAME}_alt_1.fastq.gz ${PROJ_NAME}_unpaired_1.fastq.gz ${PROJ_NAME}_alt_2.fastq.gz ${PROJ_NAME}_unpaired_2.fastq.gz ILLUMINACLIP:$PATH_TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:6:true HEADCROP:$READS1_TRIMF SLIDINGWINDOW:$QUAL_WINDOW_1:$QUAL_WINDOW_QUAL_1 MINLEN:$MIN_LENGTH_1

else
  echo "Run Trimmomatic with readtrim."
  java -jar $PATH_TRIMMOMATIC/trimmomatic-0.32.jar PE -threads 2 $READS1 $READS2 ${PROJ_NAME}_alt_1.fastq.gz ${PROJ_NAME}_unpaired_1.fastq.gz ${PROJ_NAME}_alt_2.fastq.gz ${PROJ_NAME}_unpaired_2.fastq.gz ILLUMINACLIP:$PATH_TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:6:true CROP:$TAILCROP HEADCROP:$READS1_TRIMF SLIDINGWINDOW:$QUAL_WINDOW_1:$QUAL_WINDOW_QUAL_1 MINLEN:$MIN_LENGTH_1
fi

echo "Trimmomatic finished" >>$PATH_PROJECTS_DATA/$PROJ_NAME/log.txt

# print filter infos to info file
cat run_trimmomatic.sh.e* | grep -E "Input Read Pairs" | sed 's/ Forward Only.*//' | sed 's/ *//g' | sed 's/\:/ /g' | sed 's/Both/\nBoth/g' | sed 's/(/ (/' >>$PATH_PROJECTS_DATA/$PROJ_NAME/info.txt

mkdir -p $OUTDIR_QC
$PATH_FASTQC/fastqc -o $OUTDIR_QC/ ${PROJ_NAME}_alt_1.fastq.gz
$PATH_FASTQC/fastqc -o $OUTDIR_QC/ ${PROJ_NAME}_alt_2.fastq.gz

echo "Adapter and filter steps finished. - You may proceed with mapping."
echo "FastQC finished" >>$PATH_PROJECTS_DATA/$PROJ_NAME/log.txt
