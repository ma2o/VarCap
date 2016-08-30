#!/bin/bash

# generate/reset file
>varscan.all
>lofreq1.all
>lofreq2.all
>varscan_indel_small.all
>pindel_indel_small.all
>pindel_indel_large.all
>breakdancer_indel_large.all
>delly_indel_large.all
>cortex_indel_large.all
>pindel_sv.all
>pindel_sv_del.all
>breakdancer_sv.all
>breakdancer_sv_del.all
>delly_sv.all
>delly_sv_del.all
>cortex_sv.all
>pindel_inv.all
>delly_inv.all
>cortex_inv.all

REF=$1
TOT_SNP=$( cat $REF | grep -v '#' | grep -e 'SNP' | wc -l )
TOT_INDEL_SMALL=$( cat $REF | grep -v '#' | grep -E 'DEL|INS' | grep -Ev 'DUP|ITX' | grep -E 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_INDEL_LARGE=$( cat $REF | grep -v '#' | grep -E 'DEL|INS' | grep -Ev 'DUP|ITX' | grep -Ev 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_INDEL_LARGE_INS=$( cat $REF | grep -v '#' | grep -E 'INS' | grep -Ev 'DUP|ITX' | grep -Ev 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_SV=$( cat $REF | grep -v '#' | grep -E 'DUP|ITX' | grep -Ev 'DEL' | wc -l )
TOT_SV_DEL=$( cat $REF | grep -v '#' | grep -E 'DUP|ITX' | grep -E 'DEL' | wc -l )
TOT_INV=$( cat $REF | grep -v '#' | grep -E 'INV' | wc -l )

# compute iteration number
MAX_ITER=$( ls | grep -E 'true_pos.vcf' | wc -l)

# get SNP
COUNTER=0
ITER=0
for var in $( ls | grep -E 'true_pos.vcf'); do
  if (($ITER < $MAX_ITER)); then
    cat $var | grep -e 'SNP' | grep -e 'varscan'| cut -f1 | sort -u >>varscan.all;
    cat $var | grep -e 'SNP' | grep -e 'lofreq'| grep -v 'lofreq2' | cut -f1 | sort -u >>lofreq1.all;
    cat $var | grep -e 'SNP' | grep -e 'lofreq2'| cut -f1 | sort -u >>lofreq2.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'varscan' | cut -f1 | sort -u >>varscan_indel_small.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 <= 10) print }' | cut -f1 | sort -u >>pindel_indel_small.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 > 10) print }' | cut -f1 | sort -u >>pindel_indel_large.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'breakdancer' | cut -f1 | sort -u >>breakdancer_indel_large.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'delly' | cut -f1 | sort -u >>delly_indel_large.all;
    cat $var | grep -E 'DEL|INS|IND|LI|COMPLEX' | grep -Ev 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_indel_large.all;
    
    cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>pindel_sv.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>pindel_sv_del.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>breakdancer_sv.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>breakdancer_sv_del.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>delly_sv.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>delly_sv_del.all;
    cat $var | grep -E 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_sv.all;
    
    cat $var | grep -E 'INV' | grep -e 'pindel' | cut -f1 | sort -u >>pindel_inv.all;
    cat $var | grep -E 'INV' | grep -e 'delly' | cut -f1 | sort -u >>delly_inv.all;
    cat $var | grep -E 'INV' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_inv.all;
    
    let ITER+=1
  fi
  let COUNTER+=1
done

# get counts per caller
SNP_VARSCAN=$( cat varscan.all | sort -uh | wc -l )
SNP_LOFREQ1=$( cat lofreq1.all | sort -uh | wc -l )
SNP_LOFREQ2=$( cat lofreq2.all | sort -uh | wc -l )
INDEL_VARSCAN_SMALL=$( cat varscan_indel_small.all | sort -uh | wc -l )
INDEL_PINDEL_SMALL=$( cat pindel_indel_small.all | sort -uh | wc -l )
INDEL_PINDEL_LARGE=$( cat pindel_indel_large.all | sort -uh | wc -l )
INDEL_BREAKD_LARGE=$( cat breakdancer_indel_large.all | sort -uh | wc -l )
INDEL_DELLY_LARGE=$( cat delly_indel_large.all | sort -uh | wc -l )
INDEL_CORTEX_LARGE=$( cat cortex_indel_large.all | sort -uh | wc -l )
SV_PINDEL=$( cat pindel_sv.all | sort -uh | wc -l )
SV_DEL_PINDEL=$( cat pindel_sv_del.all | sort -uh | wc -l )
SV_BREAKD=$( cat breakdancer_sv.all | sort -uh | wc -l )
SV_DEL_BREAKD=$( cat breakdancer_sv_del.all | sort -uh | wc -l )
SV_DELLY=$( cat delly_sv.all | sort -uh | wc -l )
SV_DEL_DELLY=$( cat delly_sv_del.all | sort -uh | wc -l )
SV_CORTEX=$( cat cortex_sv.all | sort -uh | wc -l )
INV_PINDEL=$( cat pindel_inv.all | sort -uh | wc -l )
INV_DELLY=$( cat delly_inv.all | sort -uh | wc -l )
INV_CORTEX=$( cat cortex_inv.all | sort -uh | wc -l )


# print variants numbers
echo "summarizing variations over iterations: $ITER of $COUNTER"
echo -e "SNP_VARSCAN\t$SNP_VARSCAN"
echo -e "SNP_LOFREQ1\t$SNP_LOFREQ1"
echo -e "SNP_LOFREQ2\t$SNP_LOFREQ2"
echo -e "INDEL_VARSCAN_SMALL\t$INDEL_VARSCAN_SMALL"
echo -e "INDEL_PINDEL_SMALL\t$INDEL_PINDEL_SMALL"
echo -e "INDEL_PINDEL_LARGE\t$INDEL_PINDEL_LARGE"
echo -e "INDEL_BREAKD_LARGE\t$INDEL_BREAKD_LARGE"
echo -e "INDEL_DELLY_LARGE\t$INDEL_DELLY_LARGE"
echo -e "INDEL_CORTEX_LARGE\t$INDEL_CORTEX_LARGE"
echo -e "SV_PINDEL\t$SV_PINDEL"
echo -e "SV_DEL_PINDEL\t$SV_DEL_PINDEL"
echo -e "SV_BREAKD\t$SV_BREAKD"
echo -e "SV_DEL_BREAKD\t$SV_DEL_BREAKD"
echo -e "SV_DELLY\t$SV_DELLY"
echo -e "SV_DEL_DELLY\t$SV_DEL_DELLY"
echo -e "SV_CORTEX\t$SV_CORTEX"
echo -e "INV_PINDEL\t$INV_PINDEL"
echo -e "INV_DELLY\t$INV_DELLY"
echo -e "INV_CORTEX\t$INV_CORTEX"

# calculate correction for multiple iterations
SNP_VARSCAN_PC=$(( ($( cat varscan.all | wc -l ) * 100) / ($TOT_SNP * $MAX_ITER) ))
SNP_LOFREQ1_PC=$(( ($( cat lofreq1.all | wc -l ) * 100) / ($TOT_SNP * $MAX_ITER) ))
SNP_LOFREQ2_PC=$(( ($( cat lofreq2.all | wc -l ) * 100) / ($TOT_SNP * $MAX_ITER) ))
INDEL_VARSCAN_SMALL_PC=$(( ($( cat varscan_indel_small.all | wc -l ) * 100) / ($TOT_INDEL_SMALL * $MAX_ITER) ))
INDEL_PINDEL_SMALL_PC=$(( ($( cat pindel_indel_small.all | wc -l ) * 100) / ($TOT_INDEL_SMALL * $MAX_ITER) ))
INDEL_PINDEL_LARGE_PC=$(( ($( cat pindel_indel_large.all | wc -l ) * 100) / ($TOT_INDEL_LARGE * $MAX_ITER) ))
INDEL_BREAKD_LARGE_PC=$(( ($( cat breakdancer_indel_large.all | wc -l ) * 100) / ($TOT_INDEL_LARGE * $MAX_ITER) ))
INDEL_DELLY_LARGE_PC=$(( ($( cat delly_indel_large.all | wc -l ) * 100) / ($TOT_INDEL_LARGE * $MAX_ITER) ))
INDEL_CORTEX_LARGE_PC=$(( ($( cat cortex_indel_large.all | wc -l ) * 100) / ($TOT_INDEL_LARGE * $MAX_ITER) ))
INDEL_CORTEX_LARGE_INS_PC=$(( ($( cat cortex_indel_large.all | wc -l ) * 100) / ($TOT_INDEL_LARGE_INS * $MAX_ITER) ))
SV_PINDEL_SV_PC=$(( ($( cat pindel_sv.all | wc -l ) * 100) / ($TOT_SV * $MAX_ITER) ))
SV_BREAKD_SV_PC=$(( ($( cat breakdancer_sv.all | wc -l ) * 100) / ($TOT_SV * $MAX_ITER) ))
SV_DELLY_SV_PC=$(( ($( cat delly_sv.all | wc -l ) * 100) / ($TOT_SV * $MAX_ITER) ))
INV_PINDEL_PC=$(( ($( cat pindel_inv.all | wc -l ) * 100) / ($TOT_INV * $MAX_ITER) ))

echo -e "Corrected rates $TOT_SNP $MAX_ITER $(cat varscan.all | wc -l)"
echo -e "SNP_VARSCAN_PC\t$SNP_VARSCAN_PC"
echo -e "SNP_LOFREQ1_PC\t$SNP_LOFREQ1_PC"
echo -e "SNP_LOFREQ2_PC\t$SNP_LOFREQ2_PC"
echo -e "INDEL_VARSCAN_SMALL_PC\t$INDEL_VARSCAN_SMALL_PC"
echo -e "INDEL_PINDEL_SMALL_PC\t$INDEL_PINDEL_SMALL_PC"
echo -e "INDEL_PINDEL_LARGE_PC\t$INDEL_PINDEL_LARGE_PC"
echo -e "INDEL_BREAKD_LARGE_PC\t$INDEL_BREAKD_LARGE_PC"
echo -e "INDEL_DELLY_LARGE_PC\t$INDEL_DELLY_LARGE_PC"
echo -e "INDEL_CORTEX_LARGE_PC\t$INDEL_CORTEX_LARGE_PC"
echo -e "INDEL_CORTEX_LARGE_INS_PC\t$INDEL_CORTEX_LARGE_INS_PC"
echo -e "SV_PINDEL_SV_PC\t$SV_PINDEL_SV_PC"
echo -e "SV_BREAKD_SV_PC\t$SV_BREAKD_SV_PC"
echo -e "SV_DELLY_SV_PC\t$SV_DELLY_SV_PC"
echo -e "INV_PINDEL_PC\t$INV_PINDEL_PC"
# echo -e "\t$"
