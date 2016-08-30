#!/bin/bash

REGEX=$1

for file in $( ls | grep -P "$REGEX" ); do
  if [ -d "${file}" ]; then
    less $file/vcfs/*filter_2.vcf | grep -e SNP | grep -Ev 'samtools|gatk|cortex' | cut -f2,10 | sort -k1,1 -un | sed 's/\:[0-9,a-z]*\:[0-9]*\:/\t/' >${file}_pct_hist.txt
    sed -i '1ipos\tcov\tfreq' ${file}_pct_hist.txt
    # generate file for rscript
    echo "mydata = read.table(\""${file}_pct_hist.txt"\", sep=\"\t\", header=TRUE, dec=\".\")" >pct_hist.R
    # echo "str(mydata)" >>pct_hist.R
    # echo "mydata <- as.numeric(mydata)" >>pct_hist.R
    echo "hist(mydata\$freq,20,xlab=\"Variant abundance\",col=\"lightgreen\",main=\"Frequencies of variant abundances\")" >>pct_hist.R
    Rscript pct_hist.R
    mv Rplots.pdf ${file}_Rplots.pdf
  fi
done



