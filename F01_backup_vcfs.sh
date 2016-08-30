#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

DESTDIR=$(echo $1 | sed 's/\/$//')

if [ "$#" -ne 1 ];then
  echo "Usage: bash E_backup_vcfs.sh <destination_directory>"
  echo ""
  exit
fi

# backup
mkdir -p ${DESTDIR}/$PROJ_NAME
cp $PATH_PROJECTS_DATA/$PROJ_NAME/* ${DESTDIR}/$PROJ_NAME/
cp -r $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs ${DESTDIR}/$PROJ_NAME/
cp -r $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw ${DESTDIR}/$PROJ_NAME/
cp -r $PATH_PROJECTS_DATA/$PROJ_NAME/annotator ${DESTDIR}/$PROJ_NAME/

# loop through files
# for file in $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf; do
# done


