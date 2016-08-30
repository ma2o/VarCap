#!/bin/bash

REGEX=$1

for file in $( ls | grep -P "$REGEX" ); do
  if [ -d "${file}" ]; then
    less $file/vcfs/*filter_2.vcf | grep -e SNP | grep -Ev 'samtools|gatk|cortex' | cut -f1,2,10 | sort -k2,2 -k1,1 -un | sed 's/\:[0-9,a-z]*\:[0-9]*\:/\t/' >${file}_pct_hist.txt
    sed -i '1ichrom\tpos\tcov\tfreq' ${file}_pct_hist.txt
    # generate file for rscript
    echo "mydata = read.table(\""${file}_pct_hist.txt"\", sep=\"\t\", header=TRUE, dec=\".\")" >pct_hist.R
    echo "mysum <- unique(mydata\$chrom)" >>pct_hist.R
    echo "attach(mtcars)" >>pct_hist.R
    echo "par(mfrow=c(2,2))" >>pct_hist.R
    echo "for (val in mysum) {sub1<-subset( mydata, chrom==val); hist(sub1\$freq,20,main=val,col=\"lightgreen\",xlab=\"Variant abundance\")}" >>pct_hist.R
    
    # echo "mydata <- as.numeric(mydata)" >>pct_hist.R
    # echo "hist(mydata\$freq,20,xlab=\"Variant abundance\",col=\"lightgreen\",main=\"Frequencies of variant abundances\")" >>pct_hist.R
    Rscript pct_hist.R
    mv Rplots.pdf ${file}_Rplots.pdf
  fi
done

# mylist <- read.table("E25_coinf_pct_hist.txt", sep="\t",header=TRUE, dec="." )
# view(mylist)
# mysum <- unique(mylist$chrom)
# subset( mylist, chrom==mysum )
# for (val in mysum) {sub1<-subset( mylist, chrom==val); print(sub1);}
# for (val in mysum) {sub1<-subset( mylist, chrom==val); hist(sub1$freq,20,main=val,col="lightgreen",xlab="Variant abundance")}
# attach(mtcars)
# par(mfrow=c(2,2))

