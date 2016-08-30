#!/bin/bash
#$-q all.q@cube[ab]*

REF_FILE=$1
VCF_FILE=$2
SNPEFF_DATA=$3
SNPEFF_PATH=$4

VCF_BASE=$(basename $VCF_FILE | sed 's/\.vcf$//')

# run snpeff (options: -no-downstream -no-upstream)
# cd $SNPEFF_PATH
# java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -i vcf -o vcf -ud 10 -c $SNPEFF_PATH/snpEff.config $SNPEFF_REF $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf
# cd $SNPEFF_PATH
java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -c $SNPEFF_PATH/snpEff.config -i vcf -o vcf -ud 0 -v $REF_FILE $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf

