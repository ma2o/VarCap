#!/bin/bash

for file in $( ls | grep -E '^[0-9]'); do less $file | sed -n '1~10p' | awk '{ print $2"\t"$3}' >${file}_10.txt; done

