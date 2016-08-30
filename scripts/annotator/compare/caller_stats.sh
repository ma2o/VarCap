#!/bin/bash

# generate/reset file
>varscan.all
>lofreq.all
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
