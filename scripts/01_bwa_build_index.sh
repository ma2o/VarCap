#!/bin/bash
#$-q all.q@cube[ab]*
#$-o mapper/bwa/bwa_indexo.log
#$-e mapper/bwa/bwa_indexe.log
#$-pe parallel 1
#$-l vf=4G

. /etc/profile

REF=$1
REFNAME=$(basename $REF | sed 's/\.f.*$//')
PATH_BWA_2=$2
PATH_INDEX=$3

mkdir bwa_index_$REFNAME
#first step: index database
${PATH_BWA_2}/bwa index -p ${PATH_INDEX}/bwa_index_$REFNAME/$REFNAME $REF

