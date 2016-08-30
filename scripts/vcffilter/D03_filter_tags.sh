#!/bin/bash

VCF=$1
OUTNAME=$( echo ${VCF%%.vcf} )
FILTER=$2
# exceptions that are also reported even if there is only support of one
EXCEPT="BP|LI|COMPLEX"

cat $VCF | awk -v regex="$FILTER" -v except="$EXCEPT" '{ if ( $7 !~ regex || $8 ~ except ) {print} }'


