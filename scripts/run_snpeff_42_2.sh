#!/bin/bash
#$-q all.q@cube[ab]*

REF_FILE=$1
VCF_FILE=$2
SNPEFF_DATA=$3
SNPEFF_PATH=$4
SNPEFFDIR=/scratch/zojer/projects/snpeff/4.2
REF_NAME=NC_005861_PAC_03

### 1. Copy and modify the config file
cp $SNPEFF_PATH/snpEff.config $SNPEFF_DATA/
sed 's|data.dir = \.\/data\/|data.dir = '"$SNPEFFDIR"'|' $SNPEFF_DATA/snpEff.config
# sed 's|data.dir = \.\/data\/|data.dir = '"$SNPEFFDIR"'|' /apps/snpeff/4.2/snpEff.config

### 2. Create SNPeff db if not existant
# 2.1 add entry to snpEff.config
HEADER=
ENTRY=
# Protochlamydia genome
$REF_NAME.genome : $REF_NAME
# add genome entry to snpEff.config
echo -e "# Protochlamydia genome" >>$SNPEFF_DATA/snpEff.config
echo -e "$REF_NAME.genome : $REF_NAME" >>$SNPEFF_DATA/snpEff.config

# 2.2 create snpeff db
cd $SNPEFFDIR
mkdir $REF_NAME
cd $REF_NAME
cp Mus_musculus.NCBIM37.61.gtf.gz genes.gtf.gz
cp Mus_musculus.NCBIM37.61.dna.toplevel.fa.gz sequences.fa

cd /path/to/snpEff
java -jar snpEff.jar build -gtf22 -v mm37.61
java -jar snpEff.jar build -gff3 -v dm5.31
java -jar /apps/snpeff/4.2/snpEff.jar build -c snpEff.config -gff3 -v NC_005861_PAC_03
java -jar /apps/snpeff/4.2/snpEff.jar build -genbank -v NC_005861_PAC_03 -c snpEff.config

### 2.a Search database
java -jar /apps/snpeff/4.2/snpEff.jar databases -c snpEff.config
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

