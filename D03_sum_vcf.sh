#!/bin/bash

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

### 6. Generate statistics
### 6.1 print coverage
# get genome coverage
# coverages are stored in vcfs_raw/*1_cov_total.txt

# get coverage avg and stdev of variants
# get SNP coverage
cd $PATH_PROJECTS_DATA/$PROJ_NAME
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*.vcf | \
grep -e 'SNP' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_snp_2.txt
# get Indel/SV coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*.vcf | \
grep -E 'DEL|INS|IND' | grep -E 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_indelsmall_2.txt
# get SV coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*.vcf | \
grep -E 'DEL|INS|INV|DUP|ITX|CTX|LI|COMPLEX' | grep -Ev 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_sv_2.txt
# get BP coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*.vcf | \
grep -E 'SVTYPE=BP' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_bp_2.txt

# print variant cov plot
# echo "library(ggplot2)" >cov_var.R
echo "totcov <- read.table(\"vcfs_raw/"${BAM_NAME_BASE}"_1_cov_total.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >cov_var.R
if [ $( cat vcfs_raw/cov_snp_2.txt | wc -l ) -gt 0 ]; then
  echo "snp <- read.table(\"vcfs_raw/cov_snp_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "snp <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat vcfs_raw/cov_indelsmall_2.txt | wc -l ) -gt 0 ]; then
  echo "indel <- read.table(\"vcfs_raw/cov_indelsmall_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "indel <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat vcfs_raw/cov_sv_2.txt | wc -l ) -gt 0 ]; then
  echo "sv <- read.table(\"vcfs_raw/cov_sv_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "sv <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat vcfs_raw/cov_bp_2.txt | wc -l ) -gt 0 ]; then
  echo "bp <- read.table(\"vcfs_raw/cov_bp_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "bp <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
echo "fname <- \""${BAM_NAME_BASE}"\"" >>cov_var.R
echo "is.list(snp)" >>cov_var.R
echo "is.list(totcov)" >>cov_var.R
echo "mysum <- unique(snp[[1]])" >>cov_var.R
echo "indelsum <- unique(indel[[1]])" >>cov_var.R
echo "totsum <- unique(totcov[[1]])" >>cov_var.R
echo "svsum <- unique(sv[[1]])" >>cov_var.R
echo "bpsum <- unique(bp[[1]])" >>cov_var.R
echo "sumlen <- length(totsum)" >>cov_var.R
echo "par(mfrow=c(3,2))" >>cov_var.R
echo "print (sumlen)" >>cov_var.R
echo "for (val in totsum) {" >>cov_var.R
echo "subtot<-subset( totcov, totcov[[1]]==val ); sub1<-subset( snp, snp[[1]]==val); subindel<-subset( indel, indel[[1]]==val);" >>cov_var.R
echo "subsv<-subset( sv, sv[[1]]==val ); subbp<-subset( bp, bp[[1]]==val);" >>cov_var.R
# plot coverage
echo "plot(subtot[c(2,3)],col=\"lightblue\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.3);" >>cov_var.R
echo "if( val %in% mysum ){ points(sub1[c(2,4)],col=\"red\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.6); mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% indelsum ){ points(subindel[c(2,4)],col=\"blue\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100));mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% svsum ){ points(subsv[c(2,4)],col=\"orange\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% bpsum ){ points(subbp[c(2,4)],col=\"black\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "par(xpd=TRUE)" >>cov_var.R
echo "legend(\"bottomleft\", inset=c(-0.18,-0.6), c(\"SNP\",\"InDel\",\"SV\",\"BP\"),lty=c(NA,NA,NA,NA),lwd=c(2.5,2.5,2.5,2.5),col=c(\"red\",\"blue\",\"orange\",\"black\"),pch=c(1,1,1,1),cex=.6)" >>cov_var.R
echo "par(xpd=FALSE)" >>cov_var.R


# plot frequency
echo "if( val %in% mysum ){ plot(sub1[c(2,5)],col=\"red\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100));mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% indelsum ){ points(subindel[c(2,5)],col=\"blue\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% svsum ){ points(subsv[c(2,5)],col=\"orange\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% bpsum ){ points(subbp[c(2,5)],col=\"black\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "}" >>cov_var.R

Rscript cov_var.R
mv Rplots.pdf vcfs/"${BAM_NAME_BASE}"_coverages.pdf

### 6.2 plot variant frequencies
# first extract frequencies, then add header
less vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*.vcf | \
grep -e SNP | grep -Ev '^#|REP|CV1' | grep -Ev 'samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[3]"\t"pct[4] }' | sort -un -k2,2 \
>vcfs_raw/${BAM_NAME_BASE}_pct_hist.txt
sed -i '1ichrom\tpos\tcov\tfreq' vcfs_raw/${BAM_NAME_BASE}_pct_hist.txt

# generate file for rscript
    echo "snp <- read.table(\""vcfs_raw/${BAM_NAME_BASE}_pct_hist.txt"\", sep=\"\t\", header=TRUE, dec=\".\")" >pct_hist.R
    echo "mysum <- unique(snp\$chrom)" >>pct_hist.R
    echo "sumlen <- length(mysum)" >>pct_hist.R
    echo "fname <- \""${BAM_NAME_BASE}"\"" >>pct_hist.R
    echo "xs = 1" >>pct_hist.R
    echo "ys = 2" >>pct_hist.R
    echo "if (sumlen > 2) { xs = 2 }" >>pct_hist.R
    echo "par(mfrow=c(ys,xs))" >>pct_hist.R
    echo "for (val in mysum) {" >>pct_hist.R
    echo "sub1<-subset( snp, chrom==val);" >>pct_hist.R
    echo "sub2 <-aggregate(cbind(cov,freq) ~ chrom+pos,data=sub1,FUN=mean, na.rm=TRUE)" >>pct_hist.R
    echo "main_name<-paste(fname, val, sep=\":\");" >>pct_hist.R
    echo "breakmod<-seq(0,100,by = 2);" >>pct_hist.R
    echo "hist(sub2\$freq[sub2\$freq<101],breaks=breakmod,main=main_name,col=\"lightgreen\",xlab=\"Variant frequency\",ylab=\"Counts\",xlim=c(0,100) )}" >>pct_hist.R

    Rscript pct_hist.R
    mv Rplots.pdf vcfs/${BAM_NAME_BASE}_frequency.pdf
