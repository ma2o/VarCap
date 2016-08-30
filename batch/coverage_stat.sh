#!/bin/bash

INDIR=$1
INNAME=$( echo ${INDIR} | sed 's/\/$//' )

# get genome coverage
echo "${INNAME}_genome_coverage"
# samtools depth $INDIR/mapper/bwa/*_v1.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "average\t",sum/NR; print "stdev\t",sqrt(sumsq/NR - (sum/NR)**2)}'

# get coverage avg and stdev of variants
echo "${INNAME}_variants_coverage"
# less $INDIR/vcfs_raw/cov_pos_1.txt | awk '{sum+=$3; sumsq+=$3*$3} END { print "average\t",sum/NR; print "stdev\t",sqrt(sumsq/NR - (sum/NR)**2)}'

# print variant cov plot
# echo "library(ggplot2)" >cov_var.R
echo "mydata = read.table(\""${INDIR}/vcfs_raw/cov_pos_1.txt"\", sep=\"\t\", header=FALSE, dec=\".\")" >cov_var.R
echo "is.list(mydata)" >>cov_var.R
echo "mysum <- unique(mydata[[1]])" >>cov_var.R
echo "sumlen = length(mysum)" >>cov_var.R
echo "attach(mtcars)" >>cov_var.R
echo "par(mfrow=c(1,1))" >>cov_var.R
echo "print (sumlen)" >>cov_var.R
echo "for (val in mysum) {sub1<-subset( mydata, mydata[[1]]==val); plot(sub1[c(2,3)],col=\"red\",main=val,xlab=\"position\",ylab=\"coverage\")}" >>cov_var.R
Rscript cov_var.R
mv Rplots.pdf "${INNAME}"_varcov_Rplots.pdf

