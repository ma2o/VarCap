#!/bin/bash

ALLIN=( "$@" )

>read_name_list.txt

for file in "${ALLIN[@]}";do
	echo -e "$file"
	TOOL=cat
	if [[ $file == *.gz ]]; then
		TOOL=zcat
	fi
	$TOOL "$file" | awk 'NR%4==1' >>read_name_list.txt
done

sort -u read_name_list.txt >read_name_list_unique.txt
