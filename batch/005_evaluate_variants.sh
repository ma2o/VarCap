#!/bin/bash

REF_VCF=$1
REGEX=$2
MAA=$3
CWD=$( pwd )
>CALLER_STATS.csv

# TN=( 14605 269163 334364 788335 1315193 1520265 1576131 1833049 1955848 2063807 2139799 )
TN=()

for file in $( ls -d */ | grep -Ev 'batch|bwa_index|reference|pdf|remap' | grep -E "$REGEX" ); do
  echo "$file MAA$MAA";
  cd $file/annotator;
  rm *vcf
  cp compare/run_compare_vcfs_stat_2.sh .
  cp compare/compare_vcfs_stat_2.pl .
  cp compare/caller_stats_5.sh .
  # calculate TP/FP
  >TP_FP.txt
  for vcf in $( ls ../vcfs | grep -e 'filter' ); do
  # add MAA filter
    NAME=$( basename ../vcfs/$vcf | sed 's/.vcf$//' )
    cat ../vcfs/$vcf | grep -v 'gatk' | awk -v maa=$MAA '{ if ($0 ~ /^#/ ) { print $0 } else { split($10,a,":"); if ( maa <= a[3] && \
		($2 != 14605 && $2 != 269163 && $2 != 334364 && $2 != 788335 && $2 != 1315193 && $2 != 1520265 && $2 != 1576131 && $2 != 1833049 && $2 != 1955848 && $2 != 2063807 && $2 != 2139799 )\
		) print } }' >MAA${MAA}_${NAME}.vcf
    bash run_compare_vcfs_stat_2.sh $REF_VCF MAA${MAA}_${NAME}.vcf >>TP_FP.txt

  done
  cat TP_FP.txt | grep -e '^[1-9]'
  # generate caller data
  bash caller_stats_5.sh $REF_VCF | grep -v '^CALLER' >>$CWD/CALLER_STATS.csv
  cd $CWD
done
# add header to collected caller stats
sed  '1i  PCT SD MRA TYPE' CALLER_STATS.csv >CALLER_STATS_H.csv

