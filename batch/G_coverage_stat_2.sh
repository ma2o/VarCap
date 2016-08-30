#!/bin/bash

INDIR=$1
BASENAME_BAM=$2

INNAME=$( echo ${INDIR} | sed 's/\/$//' )
# BASENAME_BAM=${INNAME}_all
if [ -z $BASENAME_BAM ]; then
  BASENAME_BAM=$INNAME
fi

# get genome coverage
echo "${INNAME}_genome_coverage:$INDIR/mapper/bwa/${BASENAME_BAM}_*_v1.bam"
# samtools depth $INDIR/mapper/bwa/${BASENAME_BAM}_*_v1.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "average\t",sum/NR; print "stdev\t",sqrt(sumsq/NR - (sum/NR)**2)}'
samtools depth $INDIR/mapper/bwa/${BASENAME_BAM}_*_v1.bam | awk 'NR%100==0' >$INDIR/vcfs_raw/cov_tot100_1.txt

# get coverage avg and stdev of variants
echo "${INNAME}_variants_coverage"
# less $INDIR/vcfs_raw/cov_pos_1.txt | awk '{sum+=$3; sumsq+=$3*$3} END { print "average\t",sum/NR; print "stdev\t",sqrt(sumsq/NR - (sum/NR)**2)}'
cat $INDIR/vcfs/${BASENAME_BAM}_*_1_cov_rep_sar_filter_2.vcf | grep -e 'SNP' | grep -v 'REP' | grep -Ev '^#|samtools,,,|gatk,,,|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' >$INDIR/vcfs_raw/cov_var2_1.txt
# print only multiple occurences for each position
cat $INDIR/vcfs_raw/cov_var2_1.txt | awk -F'[ \t]' '{ a[$2]++; if(a[$2] == 2) print; if (a[$2] >= 2) print }' >$INDIR/vcfs_raw/cov_var2_2.txt

# print variant cov plot
# echo "library(ggplot2)" >cov_var.R
echo "totcov <- read.table(\""${INNAME}"/vcfs_raw/cov_tot100_1.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >cov_var.R
if [ $( cat ${INNAME}/vcfs_raw/cov_var2_2.txt | wc -l ) -gt 0 ]; then
  echo "mydata <- read.table(\""${INNAME}"/vcfs_raw/cov_var2_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "mydata <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
echo "fname <- \""${INDIR}"\"" >>cov_var.R
echo "is.list(mydata)" >>cov_var.R
echo "is.list(totcov)" >>cov_var.R
echo "mysum <- unique(mydata[[1]])" >>cov_var.R
echo "totsum <- unique(totcov[[1]])" >>cov_var.R
echo "sumlen <- length(totsum)" >>cov_var.R
echo "par(mfrow=c(3,2))" >>cov_var.R
echo "print (sumlen)" >>cov_var.R
echo "for (val in totsum) {" >>cov_var.R
echo "subtot<-subset( totcov, totcov[[1]]==val ); sub1<-subset( mydata, mydata[[1]]==val);" >>cov_var.R
echo "plot(subtot[c(2,3)],col=\"lightblue\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.3);" >>cov_var.R
echo "if( val %in% mysum ){ points(sub1[c(2,4)],col=\"red\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.6) };mtext(fname,side=1,line=4,at=c(1),cex=.6);" >>cov_var.R
echo "if( val %in% mysum ){ plot(sub1[c(2,5)],col=\"red\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100));mtext(fname,side=1,line=4,at=c(1),cex=.6) }" >>cov_var.R
echo "}" >>cov_var.R

Rscript cov_var.R
mv Rplots.pdf "${BASENAME_BAM}"_varcov_Rplots.pdf

