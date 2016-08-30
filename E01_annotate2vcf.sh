#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

cd $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw

### run snpeff for variant annotation
mkdir -p $PATH_SNPEFF_DATA
echo "VARCAP:run SNPEff"

for file in $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf
do
  echo "$file"
  FILE_BASE=$( basename $file)
  FILENAME=${FILE_BASE%.vcf}
  bash $PATH_SCRIPTS/run_snpeff_42_ann.sh $SNPEFF_REF $file $PATH_SNPEFF_DATA $PATH_SNPEFF
done

