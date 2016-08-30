#!/bin/bash

# generate/reset file
>varscan.all
echo -n "varscan" >varscan.counts
>lofreq1.all
>samtools.all
echo -n "samtools" >samtools.counts
>gatk.all
echo -n "gatk" >gatk.counts
>lofreq2.all
echo -n "lofreq2" >lofreq2.counts
>varscan_indel_small.all
echo -n "varscan_indel_small" >varscan_indel_small.counts
>pindel_indel_small.all
echo -n "pindel_indel_small" >pindel_indel_small.counts
>pindel_indel_large.all
echo -n "pindel_indel_large" >pindel_indel_large.counts
>breakdancer_indel_large.all
echo -n "breakdancer_indel_large" >breakdancer_indel_large.counts
>breakdancer_indel_large_del.all
echo -n "breakdancer_indel_large_del" >breakdancer_indel_large_del.counts
>delly_indel_large.all
echo -n "delly_indel_large" >delly_indel_large.counts
>delly_indel_large_del.all
echo -n "delly_indel_large_del" >delly_indel_large_del.counts
>cortex_indel_large.all
echo -n "cortex_indel_large" >cortex_indel_large.counts
>cortex_indel_large_ins.all
echo -n "cortex_indel_large_ins" >cortex_indel_large_ins.counts
>pindel_sv.all
echo -n "pindel_sv" >pindel_sv.counts
>pindel_sv_del.all
echo -n "pindel_sv_del" >pindel_sv_del.counts
>breakdancer_sv.all
echo -n "breakdancer_sv" >breakdancer_sv.counts
>breakdancer_sv_del.all
echo -n "breakdancer_sv_del" >breakdancer_sv_del.counts
>delly_sv.all
echo -n "delly_sv" >delly_sv.counts
>delly_sv_del.all
echo -n "delly_sv_del" >delly_sv_del.counts
>cortex_sv.all
echo -n "cortex_sv" >cortex_sv.counts
>pindel_inv.all
echo -n "pindel_inv" >pindel_inv.counts
>delly_inv.all
echo -n "delly_inv" >delly_inv.counts
>cortex_inv.all
echo -n "cortex_inv" >cortex_inv.counts

REF=$1
TOT_SNP=$( cat $REF | grep -v '#' | grep -e 'SNP' | wc -l )
TOT_INDEL_SMALL=$( cat $REF | grep -v '#' | grep -E 'DEL|INS' | grep -Ev 'DUP|ITX' | grep -E 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_INDEL_LARGE=$( cat $REF | grep -v '#' | grep -E 'DEL|INS' | grep -Ev 'DUP|ITX' | grep -Ev 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_INDEL_LARGE_INS=$( cat $REF | grep -v '#' | grep -E 'INS' | grep -Ev 'DUP|ITX' | grep -Ev 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
TOT_INDEL_LARGE_DEL=$( cat $REF | grep -v '#' | grep -E 'DEL' | grep -Ev 'DUP|ITX' | grep -Ev 'LENGTH=[-]{0,1}[0-9]{1}\W' | wc -l )
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
    echo -n " $( cat $var | grep -e 'SNP' | grep -e 'varscan'| cut -f1 | sort -u | wc -l )" >>varscan.counts;
    cat $var | grep -e 'SNP' | grep -e 'samtools'| cut -f1 | sort -u >>samtools.all;
    echo -n " $( cat $var | grep -e 'SNP' | grep -e 'samtools'| cut -f1 | sort -u | wc -l )" >>samtools.counts;
    cat $var | grep -e 'SNP' | grep -e 'gatk'| cut -f1 | sort -u >>gatk.all;
    echo -n " $( cat $var | grep -e 'SNP' | grep -e 'gatk'| cut -f1 | sort -u | wc -l )" >>gatk.counts;
    cat $var | grep -e 'SNP' | grep -e 'lofreq'| grep -v 'lofreq2' | cut -f1 | sort -u >>lofreq1.all;
    cat $var | grep -e 'SNP' | grep -e 'lofreq2'| cut -f1 | sort -u >>lofreq2.all;
    echo -n " $( cat $var | grep -e 'SNP' | grep -e 'lofreq2'| cut -f1 | sort -u | wc -l )" >>lofreq2.counts;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'varscan' | cut -f1 | sort -u >>varscan_indel_small.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'varscan' | cut -f1 | sort -u | wc -l )" >>varscan_indel_small.counts;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 <= 10) print }' | cut -f1 | sort -u >>pindel_indel_small.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 <= 10) print }' | cut -f1 | sort -u | wc -l )" >>pindel_indel_small.counts;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 > 10) print }' | cut -f1 | sort -u >>pindel_indel_large.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($7 > 10) print }' | cut -f1 | sort -u | wc -l )" >>pindel_indel_large.counts;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'breakdancer' | cut -f1 | sort -u >>breakdancer_indel_large.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'breakdancer' | cut -f1 | sort -u | wc -l )" >>breakdancer_indel_large.counts;
    cat $var | grep -E 'DEL' | grep -Ev 'DUP|ITX' | grep -e 'breakdancer' | cut -f1 | sort -u >>breakdancer_indel_large_del.all;
    echo -n " $( cat $var | grep -E 'DEL' | grep -Ev 'DUP|ITX' | grep -e 'breakdancer' | cut -f1 | sort -u | wc -l )" >>breakdancer_indel_large_del.counts;
    cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'delly' | cut -f1 | sort -u >>delly_indel_large.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI' | grep -Ev 'DUP|ITX' | grep -e 'delly' | cut -f1 | sort -u | wc -l )" >>delly_indel_large.counts;
    cat $var | grep -E 'DEL' | grep -Ev 'DUP|ITX' | grep -e 'delly' | cut -f1 | sort -u >>delly_indel_large_del.all;
    echo -n " $( cat $var | grep -E 'DEL' | grep -Ev 'DUP|ITX' | grep -e 'delly' | cut -f1 | sort -u | wc -l )" >>delly_indel_large_del.counts;
    cat $var | grep -E 'DEL|INS|IND|LI|COMPLEX' | grep -Ev 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_indel_large.all;
    echo -n " $( cat $var | grep -E 'DEL|INS|IND|LI|COMPLEX' | grep -Ev 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u | wc -l )" >>cortex_indel_large.counts;
    cat $var | grep -E 'INS|COMPLEX' | grep -Ev 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_indel_large_ins.all;
    echo -n " $( cat $var | grep -E 'INS|COMPLEX' | grep -Ev 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u | wc -l )" >>cortex_indel_large_ins.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>pindel_sv.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u | wc -l )" >>pindel_sv.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>pindel_sv_del.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'pindel' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u | wc -l )" >>pindel_sv_del.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>breakdancer_sv.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u | wc -l )" >>breakdancer_sv.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>breakdancer_sv_del.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'breakdancer' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u | wc -l )" >>breakdancer_sv_del.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u >>delly_sv.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 != "DEL") print }' | cut -f1 | sort -u | wc -l )" >>delly_sv.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u >>delly_sv_del.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'delly' | awk '{ if ($6 == "DEL") print }' | cut -f1 | sort -u | wc -l )" >>delly_sv_del.counts;
    cat $var | grep -E 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_sv.all;
    echo -n " $( cat $var | grep -E 'DUP|ITX' | grep -e 'cortex' | cut -f1 | sort -u | wc -l )" >>cortex_sv.counts;
    cat $var | grep -E 'INV' | grep -e 'pindel' | cut -f1 | sort -u >>pindel_inv.all;
    echo -n " $( cat $var | grep -E 'INV' | grep -e 'pindel' | cut -f1 | sort -u | wc -l )" >>pindel_inv.counts;
    cat $var | grep -E 'INV' | grep -e 'delly' | cut -f1 | sort -u >>delly_inv.all;
    echo -n " $( cat $var | grep -E 'INV' | grep -e 'delly' | cut -f1 | sort -u | wc -l )" >>delly_inv.counts;
    cat $var | grep -E 'INV' | grep -e 'cortex' | cut -f1 | sort -u >>cortex_inv.all;
    echo -n " $( cat $var | grep -E 'INV' | grep -e 'cortex' | cut -f1 | sort -u | wc -l )" >>cortex_inv.counts;
    
    let ITER+=1
  fi
  let COUNTER+=1
done

# collect all counts to one file
awk '{ print }' *.counts >ALL.counts

# convert counts to percentages
while read -r line; do
  # echo $line
  TOOL=$( echo $line | cut -d' ' -f1 )
  TOTAL=0
  if [[ "$TOOL" == "varscan" || "$TOOL" == "lofreq2" || "$TOOL" == "samtools" || "$TOOL" == "gatk" ]];then
    TOTAL=$TOT_SNP
  elif [[ "$TOOL" == "varscan_indel_small" || "$TOOL" == "pindel_indel_small" ]]; then
    TOTAL=$TOT_INDEL_SMALL
  elif [[ "$TOOL" == "delly_indel_large" || "$TOOL" == "pindel_indel_large" || "$TOOL" == "breakdancer_indel_large" || "$TOOL" == "cortex_indel_large" ]]; then
    TOTAL=$TOT_INDEL_LARGE
  elif [[ "$TOOL" == "delly_indel_large_del" || "$TOOL" == "breakdancer_indel_large_del" ]]; then
    TOTAL=$TOT_INDEL_LARGE_DEL
  elif [[ "$TOOL" == "cortex_indel_large_ins" ]]; then
    TOTAL=$TOT_INDEL_LARGE_INS
  elif [[ "$TOOL" == "cortex_sv"  || "$TOOL" == "pindel_sv" || "$TOOL" == "delly_sv" || "$TOOL" == "breakdancer_sv" ]]; then
    TOTAL=$TOT_SV
  elif [[ "$TOOL" == "pindel_sv_del" || "$TOOL" == "delly_sv_del" || "$TOOL" == "breakdancer_sv_del" ]]; then
    TOTAL=$TOT_SV_DEL
  elif [[ "$TOOL" == "pindel_inv" || "$TOOL" == "delly_inv" || "$TOOL" == "cortex_inv" ]]; then
    TOTAL=$TOT_INV
  fi
  echo -n "$TOOL"
  echo $line | cut -d' ' -f2- | awk -v total=$TOTAL '{ s=0; for (i=1; i<=NF; i++) if ($i == 0) printf " "$i; else printf " "$i*100/total; printf "\n" }'
  
done < ALL.counts >ALL.counts.pc

# calculate average and stddev
while read -r line; do
  # echo "$line"
  TOOL=$( echo $line | cut -d' ' -f1 )
  echo -n "$TOOL "
  # calculate average
  # echo $line | cut -d' ' -f2- | awk '{ s=0; for (i=1; i<=NF; i++) s=s+$i; printf " "s/NF } END { print "\n"}'
  # echo $line | cut -d' ' -f2- | awk '{ sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf " "sum } END { print "\n"}'
  
  # average and stddev for each column
  echo $line | cut -d' ' -f2- | awk '{ sum=0; for(i=1;i<=NF;i++) { sum += $i; sumsq += ($i)^2} } END { for(i=1;i<=NR;i++) { av=sum/NF; printf "%.1f %.1f \n", av, sqrt((sumsq-sum^2/NF)/NF)} }'
done < ALL.counts.pc >ALL.counts.pc.av.std

# reformat file
# output needs to be in the folliwing format:
# caller PCT SD MRA TYPE
echo "CALLER PCT SD MRA TYPE" >ALL.counts.pc.av.std.format.csv
MRA=$( cat ../variant.config | grep -e '^MRA=' | cut -d'=' -f2 )

while read -r line; do
  TOOL=$( echo $line | cut -d' ' -f1 )
  TYPE=NONE
  if [[ $TOOL == *"_indel"* ]]; then
    TYPE=Indel
  elif [[ $TOOL == *"_sv"* || $TOOL == *"_inv"* ]]; then
    TYPE=SV
  elif [[ ! $TOOL == *"_indel"* || ! $TOOL == *"_sv"* || ! $TOOL == *"_inv"* ]]; then
    TYPE=SNP
  fi
  echo $line | awk -v mra=$MRA -v type=$TYPE '{ print $1" "$2" "$3" "mra" "type }' >>ALL.counts.pc.av.std.format.csv
done < ALL.counts.pc.av.std

cat ALL.counts.pc.av.std.format.csv
