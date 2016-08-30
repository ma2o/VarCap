#!/bin/bash

REGEX=$1
BAMPRE=$2

for file in $( ls -d */ | grep -P "$REGEX" ); do
  if [ -d "${file}" ]; then
    FILEBASE=$( echo "$file" | sed 's/\/$//' )
    BAM_BASE=${FILEBASE}_all
    # if [ -z "$BAMPRE" ]; then
    #   BAM_BASE=${file}_all
    # fi
    echo "${BAM_BASE}_bwa_100_1_cov_rep_sar_filter_2.vcf"
    less $file/vcfs/${BAM_BASE}_bwa_100_1_cov_rep_sar_filter_2.vcf | grep -e SNP | grep -Ev 'samtools,,,|gatk,,,|cortex' | cut -f1,2,10 | sort -k2,2 -k1,1 -un | sed 's/\:[0-9,a-z]*\:[0-9]*\:/\t/' >${BAM_BASE}_pct_hist.txt
    sed -i '1ichrom\tpos\tcov\tfreq' ${BAM_BASE}_pct_hist.txt
    # generate file for rscript
    echo "mydata <- read.table(\""${BAM_BASE}_pct_hist.txt"\", sep=\"\t\", header=TRUE, dec=\".\")" >pct_hist.R
    echo "mysum <- unique(mydata\$chrom)" >>pct_hist.R
    echo "sumlen <- length(mysum)" >>pct_hist.R
    echo "fname <- \""${BAM_BASE}"\"" >>pct_hist.R
    echo "xs = 1" >>pct_hist.R
    echo "ys = 2" >>pct_hist.R
    echo "if (sumlen > 2) { xs = 2 }" >>pct_hist.R
    echo "attach(mtcars)" >>pct_hist.R
    echo "par(mfrow=c(ys,xs))" >>pct_hist.R
    echo "for (val in mysum) {sub1<-subset( mydata, chrom==val); main_name<-paste(fname, val, sep=\":\"); breakmod<-seq(0,max(sub1\$freq[sub1\$freq<101]),by = 2); hist(sub1\$freq[sub1\$freq<101],breaks=breakmod,main=main_name,col=\"lightgreen\",xlab=\"Variant abundance\",xlim=c(0,100) )}" >>pct_hist.R
    
    # echo "mydata <- as.numeric(mydata)" >>pct_hist.R
    # echo "hist(mydata\$freq,50,xlab=\"Variant abundance\",col=\"lightgreen\",main=\"Frequencies of variant abundances\")" >>pct_hist.R
    Rscript pct_hist.R
    mv Rplots.pdf ${BAM_BASE}_freq_Rplots.pdf
  fi
done

# mylist <- read.table("E25_coinf_pct_hist.txt", sep="\t",header=TRUE, dec="." )
# view(mylist)
# mysum <- unique(mylist$chrom)
# sumlen <- length(mysum)
# xs = 1
# ys = 2
# if (sumlen > 2) { xs = 2 }
# subset( mylist, chrom==mysum )
# for (val in mysum) {sub1<-subset( mylist, chrom==val); print(sub1);}
# for (val in mysum) {sub1<-subset( mylist, chrom==val); hist(sub1$freq,20,main=val,col="lightgreen",xlab="Variant abundance")}
# attach(mtcars)
# par(mfrow=c(ys,xs))

