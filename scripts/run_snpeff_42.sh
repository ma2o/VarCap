#!/bin/bash
#$-q all.q@cube[ab]*

REF_FILE=$1
VCF_FILE=$2
SNPEFF_DATA=$3
SNPEFF_PATH=$4

### 1. Copy and modify the config file
cp $SNPEFF_PATH/snpEff.config $SNPEFF_DATA
sed 's/data.dir = \.\/data\//data.dir = snpeffdir/' $SNPEFF_DATA/snpEff.config

# Check for existing database
curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC005861&rettype=fasta&retmode=xml"
curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC005861&rettype=gb&retmode=text"
curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC005861&rettype=gb&retmode=text" | grep -E 'LOCUS|DEFINITION|ACCESSION|VERSION|SOURCE' >NC005861.txt
curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC005861&rettype=gb&retmode=text" | grep -E 'SOURCE' | sed -e 's/SOURCE  *//' | sed 's/ /_/g'
NC=NC005861
SOURCE=$( curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="$NC"&rettype=gb&retmode=text" | grep -E 'SOURCE' | sed -e 's/SOURCE  *//' | sed 's/ /_/g' )


### 2.a Search database
java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar databases
java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar databases | grep -e "$SOURCE"
java -Xmx4g -jar /apps/snpeff/4.2/snpEff.jar databases | grep -e "$SOURCE"

### 2. Run snpEff

VCF_BASE=$(basename $VCF_FILE | sed 's/\.vcf$//')

# run snpeff (options: -no-downstream -no-upstream)
# cd $SNPEFF_PATH
# java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -i vcf -o vcf -ud 10 -c $SNPEFF_PATH/snpEff.config $SNPEFF_REF $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf
# cd $SNPEFF_PATH
java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar eff -c $SNPEFF_DATA/snpEff.config -i vcf -o vcf -ud 0 -v $REF_FILE $VCF_FILE >$SNPEFF_DATA/$VCF_BASE.eff.vcf

