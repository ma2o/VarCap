#!/bin/bash
#$-q all.q@cube[ab]*

REF_FILE=$1
VCF_FILE=$2
SNPEFF_DATA=$3
SNPEFF_PATH=$4
SNPEFFDIR=/scratch/zojer/projects/snpeff/4.2
REF_NAME=NC_005861_PAC_03_d38_d63_d20_i46_s
REF_NAME=$REF_FILE
CH_NAME=NC_015702.1

### 1. Copy and modify the config file
cp $SNPEFF_PATH/snpEff.config $SNPEFF_DATA/
sed -i 's|data.dir = \.\/data\/|data.dir = '"$SNPEFFDIR"'|' $SNPEFF_DATA/snpEff.config
# sed 's|data.dir = \.\/data\/|data.dir = '"$SNPEFFDIR"'|' /apps/snpeff/4.2/snpEff.config

# add genome entry to snpEff.config
echo -e "# Protochlamydia genome" >>$SNPEFF_DATA/snpEff.config
echo -e "$REF_NAME.genome : $REF_NAME" >>$SNPEFF_DATA/snpEff.config
# echo -e "$CH_NAME.genome : $REF_NAME" >>$SNPEFF_DATA/snpEff.config


### 2. Run snpEff

VCF_BASE=$(basename $VCF_FILE | sed 's/\.vcf$//')

# run snpeff (options: -no-downstream -no-upstream)
# cd $SNPEFF_PATH
# java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -i vcf -o vcf -ud 10 -c $SNPEFF_PATH/snpEff.config $SNPEFF_REF $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf
# cd $SNPEFF_PATH
java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -c $SNPEFF_DATA/snpEff.config -i vcf -o vcf -ud 0 -v $REF_FILE $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf

