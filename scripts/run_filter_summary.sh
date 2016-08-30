#!/bin/bash
#$-q all.q@cube[ab]*

# usage filter_summary.sh <file>

INFILE=$1

# overview
echo "# Overview INTERGENIC/CODING, larger Indels can have both"
TOTAL_CALL=$( grep -v '#' $INFILE | wc -l )
echo -e "Total_calls\t"$TOTAL_CALL
TOTAL_LOPC=$( grep -v 'LOPC' $INFILE | grep -v '#' | wc -l )
echo -e "Total_LOPC\t"$TOTAL_LOPC
INTERGENIC=$( grep -v 'LOPC' $INFILE | grep -e 'INTERGENIC' | wc -l )
echo -e "INTERGENIC\t"$INTERGENIC
CODING=$( grep -v 'LOPC' $INFILE | grep -e 'CODING' | wc -l )
echo -e "CODING\t"$CODING
NONE=$( grep -v 'LOPC' $INFILE | grep -e 'EFF=NONE' | wc -l )
echo -e "NONE\t"$NONE

# snps
echo "# Overview SNPs"
SNPS=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | wc -l )
echo -e "SNP_TOTAL\t"$SNPS
SNPS_CODING=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -v 'INTERGENIC' | grep -e 'CODING' | wc -l )
SNPS_NONCOD=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -e 'INTERGENIC' | wc -l )
echo -e "SNP_NONCOD\t"$SNPS_NONCOD
echo -e "SNP_CODING\t"$SNPS_CODING
SYN_NONSYN=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -e 'SYNONYMOUS' | wc -l )
NONSYN=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -e 'EFF=NON_SYNONYMOUS' | wc -l )
SYN=$( grep -v 'LOPC' $INFILE | grep -e 'EFF=SYNONYMOUS' | wc -l )
echo -e "SNPS_NONSYN\t"$NONSYN
SNPS_STOP=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -e 'EFF=STOP' | wc -l )
echo -e "SNPS_STOP\t"$SNPS_STOP
SNPS_START=$( grep -v 'LOPC' $INFILE | grep -e 'SNP' | grep -e 'EFF=START' | wc -l )
echo -e "SNPS_START\t"$SNPS_START
echo -e "SNPS_SYN\t"$SYN

# indels (large Indels can be CODING and INTERGENIC)
echo "# Overview Indels - large Indels can be CODING and INTERGENIC"
INDEL=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | wc -l )
echo -e "INDEL_TOTAL\t"$INDEL
INDEL_INT=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | grep -e 'INTERGENIC' | wc -l )
echo -e "INDEL_INT\t"$INDEL_INT
INDEL_COD=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | grep -e 'CODING' | wc -l )
echo -e "INDEL_COD\t"$INDEL_COD
INDEL_FS=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | grep -e 'CODING' | grep -e 'FRAME_SHIFT' | wc -l )
echo -e "INDEL_FS\t"$INDEL_FS
INDEL_STOP=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | grep -e 'CODING' | grep -e 'STOP_' | wc -l )
echo -e "INDEL_STOP\t"$INDEL_STOP
INDEL_START=$( grep -v 'LOPC' $INFILE | grep -E 'DEL|INS|IND' | grep -e 'CODING' | grep -e 'START_' | wc -l )
echo -e "INDEL_START\t"$INDEL_START

# BRP
echo "# Overview BRP - are BP,ITX,DUP"
BRP=$( grep -v 'LOPC' $INFILE | grep -e 'BRP' | wc -l )
echo -e "BRP_TOTAL\t"$BRP

# misc cortex calls
echo "#  Misc cortex calls"
PH_S=$( grep -v 'LOPC' $INFILE | grep -e 'PH_S' | wc -l )
COMPLEX=$( grep -v 'LOPC' $INFILE | grep -e 'COMPLEX' | wc -l )
echo -e "PH_S\t"$PH_S
echo -e "COMPLEX\t"$COMPLEX

# total variants corrected
echo "# Total variants corrected"
TOTAL_VAR=$(( $SNPS + $INDEL + $BRP + $COMPLEX ))
echo -e "VARIANTS_CORR\t"$TOTAL_VAR
