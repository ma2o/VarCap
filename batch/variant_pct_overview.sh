#!/bin/bash

REGEX=$1

declare -a PCT=(0 0 0 0 0 0 0 0 0 0)

for file in $( ls | grep -P "$REGEX" ); do
  echo $file
  less $file/vcfs/*filter_2.vcf | grep -e SNP | grep -Ev 'samtools|gatk|cortex' | cut -f2,10 | sort -k1,1 -un | sed 's/\:[0-9,a-z]*\:[0-9]*\:/\t/' | wc -l
  # less $file/vcfs/*filter_2.vcf | grep -e SNP | grep -Ev 'samtools|gatk|cortex' | cut -f2,10 | sort -k1,1 -un | sed 's/\:[0-9,a-z]*\:[0-9]*\:/\t/' >${file}_pct.txt
  while read line; do
    PCT_VAL=$( echo "$line" | cut -f3 )
    PCT_VALT=${PCT_VAL%.*}
    # echo $PCT_VALT
    if [[ "$PCT_VALT" -ge 90 ]]; then
      let PCT[9]++
    elif [[ "$PCT_VALT" -ge 80 && "$PCT_VALT" -lt 90 ]]; then
      let PCT[8]++
    elif [[ "$PCT_VALT" -ge 70 && "$PCT_VALT" -lt 80 ]]; then
      let PCT[7]++
      echo $line
    elif [[ "$PCT_VALT" -ge 60 && "$PCT_VALT" -lt 70 ]]; then
      let PCT[6]++
      echo $line
    elif [[ "$PCT_VALT" -ge 50 && "$PCT_VALT" -lt 60 ]]; then
      let PCT[5]++
      echo $line
    elif [[ "$PCT_VALT" -ge 40 && "$PCT_VALT" -lt 50 ]]; then
      let PCT[4]++
      echo $line
    elif [[ "$PCT_VALT" -ge 30 && "$PCT_VALT" -lt 40 ]]; then
      let PCT[3]++
      echo $line
    elif [[ "$PCT_VALT" -ge 20 && "$PCT_VALT" -lt 30 ]]; then
      let PCT[2]++
      echo $line
    elif [[ "$PCT_VALT" -ge 10 && "$PCT_VALT" -lt 20 ]]; then
      let PCT[1]++
      echo $line
    elif [[ "$PCT_VALT" -ge 2 && "$PCT_VALT" -lt 10  ]]; then
      let PCT[0]++
      echo $line
    fi
  done <"${file}_pct.txt"
  echo ${PCT[@]}
  PCT=(0 0 0 0 0 0 0 0 0 0)
  # sed -i '1ipos\tfreq' *_all.txt
done

