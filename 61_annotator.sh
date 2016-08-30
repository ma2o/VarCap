#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

VCF_FILE2ANN=$1

### runs collect_variants_varcap.pl to collect variants from the different callers into one vcf
PATH_TO_VARIANTS_PERL=$PATH_SCRIPTS/annotator
# export PATH_CALLER
RC_RATIO=$(((100 * ${SUBSAMPLE_SIZE_ALT}) / (${SUBSAMPLE_SIZE_REF} + ${SUBSAMPLE_SIZE_ALT}) | bc))
FILENAME_BASE=${FILENAME_BASE}_${RC_RATIO}_v1
echo $FILENAME_BASE
cd $PATH_PROJECTS_DATA/$PROJ_NAME

### run snpeff for variant annotation
mkdir $PATH_SNPEFF_DATA
echo "VARCAP:run SNPEff"
sh $PATH_SCRIPTS/run_snpeff.sh $VCF_FILE2ANN $PATH_SNPEFF_DATA $PATH_SNPEFF $PATH_VCF


