#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo ""
  echo "Usage: bash 00_prepare.sh <path/source_dir> <project_path/project_name> [optional: reference]"
  echo ""
  echo "The folder <source_dir> hosts all source files with following formats: fastq,fastq.gz,bam"
  echo "The files can be present either all within the <source_dir>, or each one within their own subfolder."
  echo "The folder <project_name> will be generated within the path of <project_path>."
  echo "All necessary files will be copied to folder <project_name> and should be run from there."
  echo "Additional optional input: bash 00_prepare.sh <data_dir> <project_name> [path/ref.fasta] "
  echo ""
  exit
fi

SOURCE_DIR=$( echo "$1" | sed 's/\/$//')
PROJ_DIR=$( echo "$2" | sed 's/\/$//')
PROJ_PATH=$( echo "${PROJ_DIR%/*}" | sed 's/\/$//' )
PROJ_NAME=$( echo "${PROJ_DIR##*/}" | sed 's/\/$//' )
REF=$3

# required paths
PATH_VARCAP=$( pwd | sed 's/\/$//' )
# check if we execute in correct dir and load program variables, otherwise exit
if [[ ! -e program.config ]]; then
  echo "Wrong execution folder: $PATH_VARCAP"
  exit 1
fi
# load program variables
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# print input data/paths
echo "$PROJ_PATH : $PROJ_NAME : $SOURCE_DIR : $REF"

# create varcap data dir ( proj dir )
mkdir -p $PROJ_DIR
if [[ -d $PROJ_DIR ]]; then
    echo "Project directory:$PROJ_DIR created."
elif [[ ! -d $PROJ_DIR ]]; then
  echo "Project directory:$PROJ_DIR could not be created";
  exit 1;
fi

# run the varcapproject initiation
cd $SOURCE_DIR

# function to get project name depending on filestructure and input file format
# -either files are all in one folder, or each file in seperate folder
# -format is either fastq, fastq.gz or bam
project_name() {
  local INFILE=$1
  local INFILE_MOD=NA
  if [[ -d $INFILE ]]; then
    INFILE_MOD=$( echo $INFILE | sed 's/\/$//' )
    echo "$INFILE_MOD"
  elif [[ -f $INFILE ]]; then
      if [[ $INFILE == *bam ]]; then
        # echo "bam file"
        INFILE_MOD=$( echo $INFILE | sed -e 's/[_-\/]\{0,1\}[12]\{0,1\}\.bam$//' )
      elif [[ $INFILE == *gz ]]; then
        # echo "gzipped fastq file"
        INFILE_MOD=$( echo $INFILE | sed -e 's/[_-\/]\{0,1\}[12]\{0,1\}\.f.*q.gz$//' )
      elif [[ $INFILE == *fastq ]] || [[ $INFILE == *fq ]]; then
        # echo "fastq file"
        INFILE_MOD=$( echo $INFILE | sed -e 's/[_-\/]\{0,1\}[12]\{0,1\}\.f.*q$//' )
      fi
    echo "$INFILE_MOD"
  fi
}

for file in $( ls ); do
  # $( ls | grep -E "fq$|fastq$|fastq.gz$|bam$" | grep -Ev "2.fq$|2.fastq$|2.fastq.gz$|2.fq.gz$" )
  ### 0. get experiment name to create project/experiment folders
  echo "FILE:$file"
  cd $SOURCE_DIR
  EXP_NAME=$( project_name $file )
  # EXP_NAME=$( echo $file | sed 's/\/$//' )
  echo "$EXP_NAME"
  
  ### 1. create project/experiment folders ###
  # check if data output folder (varcap_data) exists, else create it
  mkdir -p $PROJ_PATH/$PROJ_NAME/
  mkdir -p $PROJ_PATH/$PROJ_NAME/batch
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME/filter
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME/mapper
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME/caller
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME/annotator
  mkdir -p $PROJ_PATH/${PROJ_NAME}/$EXP_NAME/logs
  
  # copy config file and top level scripts to project folder, to store run information
  echo "Copy main scripts to project folder."
  cp -u ${PATH_VARCAP}/variant.config ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}/variant.config
  cp -u ${PATH_VARCAP}/A??_*.sh ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  cp -u ${PATH_VARCAP}/B??_*.sh ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  cp -u ${PATH_VARCAP}/C??_*.sh ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  cp -u ${PATH_VARCAP}/D??_*.sh ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  cp -u ${PATH_VARCAP}/H??_*.sh ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  cp -ru ${PATH_VARCAP}/scripts/annotator/compare ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}/annotator
  cp -ru ${PATH_VARCAP}/batch/A*.sh ${PROJ_PATH}/${PROJ_NAME}/batch
  cp -ru ${PATH_VARCAP}/batch/00*.sh ${PROJ_PATH}/${PROJ_NAME}
  
  ### 2. create symbolic links to data
  # check if we are in data folder or top folder
  INFILES_PATH=NA
  # INFILES=NA
  if [[ -d $file  ]]; then
    cd $file
    INFILES_PATH=$( pwd | sed 's/\/$//')
    # create link to raw input files in /raw folder
    mkdir -p ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw
    cd ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw
    for file2 in $( find $INFILES_PATH/* | grep -E 'fq$|fastq$|gz$|bam$' | sed 's/.*\///' ); do
      # echo "$INFILES_PATH"
      # echo "$file2"
      ln -s $INFILES_PATH/$file2 ${EXP_NAME}_$file2
    done
    # INFILES=$( ls | grep -E 'fq$|fastq$|gz$|bam$' )
  elif  [[ -f $file ]]; then
    INFILES_PATH=$( pwd | sed 's/\/$//')
    # create link to raw input files in /raw folder
    mkdir -p ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw
    cd ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw
    for file2 in $( find $INFILES_PATH/$EXP_NAME* | grep -E 'fq$|fastq$|gz$|bam$' | sed 's/.*\///' ); do
      ln -s $INFILES_PATH/$file2 $file2
    done
  fi
  

  ### 3. create/update variant.config file
  # determine readlength
  READLEN_MAX="max_length\t100"
  INFILE_END=$( find ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/* | grep -E 'fq$|fastq$|gz$|bam$' | head -n1 )
  # echo "$INFILE_END"
  if [[ $INFILE_END == *bam ]]; then
    # echo "readlen bam"
    #READLEN=$( $PATH_SAMTOOLS/samtools view ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/*bam | head -n10000 | awk '{print length($10)}' | sort | uniq -c | sort -nr | head -n1 | cut -f2 )
    READLEN=$( $PATH_SAMTOOLS/samtools view ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/*bam | head -n10000 | \
      awk '{len=length($10);sum+=len;if(min==""){min=max=len}; if(len>max) {max=len}; if(len< min) {min=len};} END {print "average\t",sum/NR;print "max_length\t",max;print "min_length\t",min;}')
    READLEN_MAX=$( echo -e "$READLEN" | grep -e 'max_length' | cut -f2 | sed 's/ //g' )
  elif [[ $INFILE_END == *gz ]]; then
    # echo "readlen gz"
    #READLEN=$( zcat ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/*1.f*q.gz | grep -E -A1 '^@' | head -n10000 | grep -v @ | grep -v - | awk '{ print length($0); }' | sort | uniq -c | sort -nr | head -n1 | cut -f2 )
    READLEN=$( zcat ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/*1.f*q.gz | grep -E -A1 '^@' | head -n10000 | grep -v @ | grep -v - | \
      awk '{len=length($0);sum+=len;if(min==""){min=max=len}; if(len>max) {max=len}; if(len< min) {min=len};} END {print "average\t",sum/NR;print "max_length\t",max;print "min_length\t",min;}')
    READLEN_MAX=$( echo -e "$READLEN" | grep -e 'max_length' | cut -f2 | sed 's/ //g' )
  elif [[ $INFILE_END == *fastq ]] || [[ $INFILE_END == *fq ]]; then
    # echo "readlen fastq"
    READLEN=$( cat ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/raw/*1.f*q | grep -E -A1 '^@' | head -n10000 | grep -v @ | grep -v - | \
      awk '{len=length($0);sum+=len;if(min==""){min=max=len}; if(len>max) {max=len}; if(len< min) {min=len};} END {print "average\t",sum/NR;print "max_length\t",max;print "min_length\t",min;}')
    READLEN_MAX=$( echo -e "$READLEN" | grep -e 'max_length' | cut -f2 | sed 's/ //g' )
  fi
  
  # create info.txt file
  echo -e "name\t$EXP_NAME" >${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}/info.txt
  echo -e "$READLEN" | sed -e 's/  *//g' >>${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}/info.txt
  READ_LEN=${READLEN_MAX}
  SUBSIZEALT=-1
  REP=1
  IS_SIZE=250
  MIN_LENGTH_1=40
  MIN_LENGTH_2=40
  
  # edit variant config to add varcap path and project name
  cd ${PROJ_PATH}/${PROJ_NAME}/${EXP_NAME}
  # edit varcap path
  sed -i 's#^PATH_VARCAP=.*$#PATH_VARCAP='"${PATH_VARCAP}"'#' variant.config
  # edit project path
  sed -i 's#^PATH_PROJECTS_DATA=.*$#PATH_PROJECTS_DATA='"${PROJ_PATH}/${PROJ_NAME}"'#' variant.config
  # edit project name and bam name and bam input
  sed -i 's#^PROJ_NAME=.*$#PROJ_NAME='"${EXP_NAME}"'#' variant.config
  sed -i 's#^BAM_NAME_BASE=.*$#BAM_NAME_BASE='"${EXP_NAME}"'#' variant.config
  # edit settings
  sed -i 's#^SUBSAMPLE_SIZE_ALT=.*$#SUBSAMPLE_SIZE_ALT='"${SUBSIZEALT}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  sed -i 's#^REPEATS=.*$#REPEATS='"${REP}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  sed -i 's#^READLENGTH=.*$#READLENGTH='"${READ_LEN}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  sed -i 's#^INSERT_SIZE=.*$#INSERT_SIZE='"${IS_SIZE}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  sed -i 's#^MIN_LENGTH_1=.*$#MIN_LENGTH_1='"${MIN_LENGTH_1}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  sed -i 's#^MIN_LENGTH_2=.*$#MIN_LENGTH_2='"${MIN_LENGTH_2}"'#' ${PROJ_PATH}/$PROJ_NAME/$EXP_NAME/variant.config
  
  # process reference
  if [ ! -z "$REF" ]; then
    sed -i 's#^REF_FA_ORIGINAL=.*$#REF_FA_ORIGINAL='"${REF}"'#' variant.config;
  fi

  # report log
  echo -e "1. Project preperation finished:${PROJ_PATH}/$PROJ_NAME/$EXP_NAME" >logs/log.txt
  
done

# report end
echo -e "1. Project preperation finished. - Change to project folder (${PROJ_PATH}/$PROJ_NAME) to proceed."

