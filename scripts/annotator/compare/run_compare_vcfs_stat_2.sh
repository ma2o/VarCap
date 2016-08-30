#!/bin/bash

REF=$1
ALT=$2
# REF=/scratch/zojer/projects/variator_data/09_135Av07_Proto_am_U25.vcf
# ALT=/scratch/zojer/projects/varcap_data/E25_pirs_500x_250is_10std_09_135Av07_Proto_am_U25/vcfs/E25_pirs_500x_250is_10std_135Av07_Proto_am_U25_bwa_3_1_cov_rep_sar_filter_2.vcf

if [ "$#" -lt 2 ];then
  echo "Usage: bash run_compare_vcfs_stat_2.sh <ref.vcf> <alt.vcf> .."
  echo "The ref.vcf contains the known variants, whereas the alt.vcf is evaluated (against ref.vcf) into TP,FN,FP."
  exit
fi

# get name
REF_NAME=$( basename $REF | sed 's/.vcf$//' )
ALT_NAME=$( basename $ALT | sed 's/.vcf$//' )
OUTNAME=${REF_NAME}_${ALT_NAME}
echo $OUTNAME

# run comparison
perl compare_vcfs_stat_2.pl $REF $ALT

TOTAL=$( less $REF | grep -v '#' | cut -f2 | sort -u | wc -l )
TOT_SNP=$( less $REF | grep -v '#' | grep -e 'SNP' | cut -f2 | sort -u | wc -l )
TOT_IND=$( less $REF | grep -v '#' | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | cut -f2 | sort -u | wc -l )
TOT_DUX=$( less $REF | grep -v '#' | grep -E 'DUP|ITX' | cut -f2 | sort -u | wc -l )
TOT_INV=$( less $REF | grep -v '#' | grep -e 'INV' | cut -f2 | sort -u | wc -l )

TP=$( less ${OUTNAME}.true_pos.vcf | grep -v '#' | cut -f1 | sort -u | wc -l )
TP_SNP=$( less ${OUTNAME}.true_pos.vcf | grep -v '#' | grep -e 'SNP' | cut -f1 | sort -u | wc -l )
TP_IND=$( less $OUTNAME.true_pos.vcf | grep -v '#' | cut -f1,2,3 | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | cut -f1 | sort -u | wc -l )
TP_DUX=$( less $OUTNAME.true_pos.vcf | grep -v '#' | grep -E 'DUP|ITX' | cut -f1 | sort -u | wc -l )
TP_INV=$( less $OUTNAME.true_pos.vcf | grep -v '#' | grep -E 'INV' | cut -f1 | sort -u | wc -l )

FN=$( less $OUTNAME.false_neg.vcf | grep -v '#' | cut -f2 | sort -u | wc -l )
FN_SNP=$( less $OUTNAME.false_neg.vcf | grep -v '#' | grep -e 'SNP' | cut -f2 | sort -u | wc -l )
FN_IND=$( less $OUTNAME.false_neg.vcf | grep -v '#' | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | cut -f2 | sort -u | wc -l )
FN_DUX=$( less $OUTNAME.false_neg.vcf | grep -v '#' | grep -E 'DUP|ITX' | cut -f2 | sort -u | wc -l )
FN_INV=$( less $OUTNAME.false_neg.vcf | grep -v '#' | grep -E 'INV' | cut -f2 | sort -u | wc -l )

FP=$( less $OUTNAME.false_pos.vcf | grep -v '#' | cut -f2 | sort -u | wc -l )
FP_SNP=$( less $OUTNAME.false_pos.vcf | grep -v '#' | grep -e 'SNP' | cut -f2 | sort -u | wc -l )
FP_IND=$( less $OUTNAME.false_pos.vcf | grep -v '#' | grep -E 'DEL|INS|IND|LI' | cut -f2 | sort -u | wc -l )
FP_DUX=$( less $OUTNAME.false_pos.vcf | grep -v '#' | grep -E 'DUP|ITX' | cut -f2 | sort -u | wc -l )
FP_INV=$( less $OUTNAME.false_pos.vcf | grep -v '#' | grep -E 'INV' | cut -f2 | sort -u | wc -l )
FP_COM=$( less $OUTNAME.false_pos.vcf | grep -v '#' | grep -E 'COMPLEX' | cut -f2 | sort -u | wc -l )

# calculate percentages
PTP=$( echo "scale=2; $TP * 100 / $TOTAL" | bc )
PTP_SNP=$( echo "scale=2; $TP_SNP * 100 / $TOT_SNP" | bc )
PTP_IND=$( echo "scale=2; $TP_IND * 100 / $TOT_IND" | bc )
PTP_DUX=$( echo "scale=2; $TP_DUX * 100 / $TOT_DUX" | bc )
PTP_INV=$( echo "scale=2; $TP_INV * 100 / $TOT_INV" | bc )

echo "TOTAL TOT_SNP TOT_IND TOT_DUX TOT_INV TP PTP FN FP TP_SNP PTP_SNP TP_IND PTP_IND TP_DUX PTP_DUX TP_INV PTP_INV FN_SNP FN_IND FN_DUX FN_INV FP_SNP FP_IND FP_DUX FP_INV FP_COM"
echo "$TOTAL $TOT_SNP $TOT_IND $TOT_DUX $TOT_INV $TP $PTP $FN $FP $TP_SNP $PTP_SNP $TP_IND $PTP_IND $TP_DUX $PTP_DUX $TP_INV $PTP_INV $FN_SNP $FN_IND $FN_DUX $FN_INV $FP_SNP $FP_IND $FP_DUX $FP_INV $FP_COM"

