#!/bin/bash

REGEX=$1

for file in $( ls -d */ | grep -E "$REGEX" ); do
	FNRAW=$( echo $file | sed 's/_[0-9\/]*$//' )
	echo -en "$FNRAW\t";
	cat ../$FNRAW/reference/*summary.txt | grep -e "topref" | tail -n 1;
	echo -en "${file%/}\t"; 
	cat $file/reference/*summary.txt | grep -e "topref" | tail -n 1; 


done | cut -f1,3,8 | sort -u -k1,1 >comp_R1_R2_all.txt

# sample list
cat comp_R1_R2_all.txt | cut -f1 | sed 's/_[0-9]*$//' | sort -u >sample_list_all.txt

cat sample_list_all.txt | while read line; do
	echo -en "$line\t"
	grep -e "$line" comp_R1_R2_all.txt | cut -f2 | sed 's/uid.*$//;s/Chlamydia_trachomatis_/Ctr_/;s/_*$//' | tr "\n" "\t" | sed 's/\t$/\n/'
	
done
