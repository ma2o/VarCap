#!/bin/bash

# generate/reset file
>snp.all
>indel.all
>dux.all
>inv.all

# maximal iterations
if [ "$#" -lt 1 ];then
  echo "Usage: combine_iteration_stats.sh <iterations>"
  exit
fi
MAX_ITER=$1

# get SNP
COUNTER=0
ITER=0
for var in $( ls | grep -E 'true_pos.vcf'); do
  if (($ITER < $MAX_ITER)); then
    cat $var | grep -e 'SNP' | cut -f1 | sort -u >>snp.all;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | cut -f1 | sort -u >>indel.all;
    cat $var | grep -E 'DUP|ITX' | cut -f1 | sort -u >>dux.all;
    cat $var | grep -E 'INV' | cut -f1 | sort -u >>inv.all;
    let ITER+=1
  fi
  let COUNTER+=1
done
# get frequency of variations
SNP_ITER=$( cat snp.all | sort -n | uniq -c | sort -k1,1n -k2,2n | sed 's/^ *//g' )
SNP_ITER_FREQ=$( echo "$SNP_ITER" | sed 's/^ *//g' | cut -d' ' -f1 | uniq -c | sed 's/^ *//g' )
INDEL_ITER=$( cat indel.all | sort -n | uniq -c | sort -k1,1n -k2,2n | sed 's/^ *//g' )
INDEL_ITER_FREQ=$( echo "$INDEL_ITER" | sed 's/^ *//g' | cut -d' ' -f1 | uniq -c | sed 's/^ *//g' )
DUX_ITER=$( cat dux.all | sort -n | uniq -c | sort -k1,1n -k2,2n | sed 's/^ *//g' )
DUX_ITER_FREQ=$( echo "$DUX_ITER" | sed 's/^ *//g' | cut -d' ' -f1 | uniq -c | sed 's/^ *//g' )
INV_ITER=$( cat inv.all | sort -n | uniq -c | sort -k1,1n -k2,2n | sed 's/^ *//g' )
INV_ITER_FREQ=$( echo "$INV_ITER" | sed 's/^ *//g' | cut -d' ' -f1 | uniq -c | sed 's/^ *//g' )

# SNP_ITER=$( cat snp.all | sort -uh | wc -l )
# INDEL_ITER=$( cat indel.all | sort -uh | wc -l )
# DUX_ITER=$( cat dux.all | sort -uh | wc -l )
# INV_ITER=$( cat inv.all | sort -uh | wc -l )

# write position frequency to file
echo "$SNP_ITER">snp_pos_freq.txt
echo "$INDEL_ITER">indel_pos_freq.txt
echo "$DUX_ITER">dux_pos_freq.txt
echo "$INV_ITER">inv_pos_freq.txt

# print variants numbers
echo "summarizing variations over iterations: $ITER of $COUNTER"
echo "print the frequency of observed variants: variant_count frequency"
echo "SNP"
echo "$SNP_ITER_FREQ"
echo "INDEL"
echo "$INDEL_ITER_FREQ"
echo "DUX"
echo "$DUX_ITER_FREQ"
echo "INV"
echo "$INV_ITER_FREQ"

