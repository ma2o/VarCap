#!/bin/bash

>mapping_info.txt
# for file in $( ls -d */ | grep -e 'ERR'); do echo $file; samtools flagstat $file/mapper/bwa/*.bam | head -n 3 | grep -v 'duplicates' >>mapping_info.txt; done
for file in $( ls -d */ | grep -e 'ERR'); do echo $file >>mapping_info.txt; samtools flagstat $file/mapper/bwa/*all_bwa_100_v1.bam | grep -E 'total|mapped' >>mapping_info.txt; done
cat mapping_info.txt | grep -A2 'ERR' | grep -Eo "^ERR.*|^[0-9]*|mapped (.*)" | tr "\n" "\t" | sed 's/)\t/)\n/g' >mapping_info_sort.txt
for file in $( ls -d */ | grep -E 'ERR'); do echo -n "$file "; cat $file/info.txt | sed -e 's/  */ /g' | tail -n3 | cut -d' ' -f2,3 | tr "\n" "\t" | sed 's/)\t/)\n/g' | sed 's/ /\t/g'; done >filter_info.txt

echo -e "RL      TOTCOV  FILCOV  FILPCT  MAPTOT  MAPMAP  IN      MAPPCT" >filter_mapped_info_sort.txt
join <( sort -k1,1 filter_info.txt) <( sort -k1,1 mapping_info_sort.txt) | sed 's/[(,),%]//g' | sed 's/:-nan//g' >>filter_mapped_info_sort.txt
cat filter_mapped_info_sort.txt | sed -e 's/  */ /g' | tr " " "\t" >filter_mapped_info_sort_tab.txt
mv filter_mapped_info_sort_tab.txt filter_mapped_info_sort.txt
