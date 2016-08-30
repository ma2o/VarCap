#!/bin/bash

BAMFILE=$1

samtools view -f66 $BAMFILE | cut -f 9 | sed 's/^-//' >InsertSizeMetrics_samtools.txt
grep -v '[0-9]\{5,\}' InsertSizeMetrics_samtools.txt >InsertSizeMetrics_samtools_1k.txt
Rscript insert_size.R >>IS_sum.txt

