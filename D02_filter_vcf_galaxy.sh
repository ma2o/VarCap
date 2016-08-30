#!/bin/bash
#
#
# Main script for Varcap filter
#
#

PATH_CALLER_COLLECT=$1
REF_FA=$2
PATH_BAM=$3
BAM_NAME_BASE=varcap
MIN_CPV=1
REPEATS=1
PATH_SAMTOOLS='/opt/apps/varcap_filter/2.7.1_galaxy/samtools/0.1.19'
PATH_VMATCH='/opt/apps/varcap_filter/2.7.1_galaxy/vmatch/2.2.4'
PATH_TO_VARIANTS_PERL='/opt/apps/varcap_filter/2.7.1_galaxy/vcffilter'
INSERT_SIZE=$4
MRA=$5
MAA=$6

###
###
### Get raw output files of all callers and join them into a master VCF
###
###
perl -I $PATH_TO_VARIANTS_PERL $PATH_TO_VARIANTS_PERL/collect_variants_varcap_2vcf_02.pl $REF_FA $REPEATS $PATH_CALLER_COLLECT $BAM_NAME_BASE

# convert cortex chrom to ref chrom name

  #>${BAM_NAME_BASE}.vcf.mod

  #CHROMS=$( cat "${BAM_NAME_BASE}.vcf" | grep -Ev '^#|cortex' | cut -f1 | sort -u )
  #cat "${BAM_NAME_BASE}.vcf" | while read -r line; do
    #if [[ $line == *cortex* ]]; then
      #CCHR=$( echo -e "$line" | cut -f1 )
      #CREST=$( echo -e "$line" | cut -f2- )
      #NCHR=$( echo -e "$CHROMS" | grep -e "$CCHR" )
      #echo -e "$NCHR\t$CREST" >>${BAM_NAME_BASE}.vcf.mod
    #else
      #echo -e "$line" >>${BAM_NAME_BASE}.vcf.mod
    #fi
  #done
   #sort according to chrom and pos
  #>HEAD.txt
  #>BODY.txt
  #cat ${BAM_NAME_BASE}.vcf.mod | while read -r line; do
    #if [[ $line ==  \#* ]]; then
      #echo -e "$line" >>HEAD.txt
    #else
      #echo -e "$line" >>BODY.txt
    #fi
  #done
  #cat BODY.txt | sort -k1,1 -k2,2n >BODY.sort.txt
  #cat HEAD.txt BODY.sort.txt >${BAM_NAME_BASE}.vcf.mod
  #rm HEAD.txt BODY.sort.txt BODY.txt
  #mv ${BAM_NAME_BASE}.vcf.mod ${BAM_NAME_BASE}.vcf

###
###
### Calculate coverages and percentages for positions and add them to vcf
###
###
$PATH_SAMTOOLS/samtools depth $PATH_BAM | awk 'NR%500==0' >${BAM_NAME_BASE}_cov_total.txt
cat ${BAM_NAME_BASE}_cov_total.txt | awk '{ a[$1]+=$3; b[$1]=NR; next } END {for(i in a) { av_cov=(a[i]/b[i]); print i"\t"av_cov } }' >${BAM_NAME_BASE}_cov_av.txt
grep -v '#' ${BAM_NAME_BASE}*.vcf | cut -f 1,2 >pos.bed
$PATH_SAMTOOLS/samtools depth -b pos.bed ${PATH_BAM} >${BAM_NAME_BASE}_cov_pos.txt
perl -I $PATH_TO_VARIANTS_PERL $PATH_TO_VARIANTS_PERL/get_coverage_2vcf2.pl "${BAM_NAME_BASE}_1.vcf" "${BAM_NAME_BASE}_cov_pos.txt" "${BAM_NAME_BASE}_cov_av.txt" >${BAM_NAME_BASE}_cov.vcf

###
###
### Search for repetitive elements within the reference genome that are longer than insert size and tag homopolymers
###
###
echo "VARCAP: search for repetitive elements using vmatch"
mkdir -p vmatch/mkvtree
cd vmatch/mkvtree
REF_IDX_NAME=$( basename $REF_FA )
OUT_NAME_BASE=$( echo $REF_IDX_NAME | sed 's/\..*$//' )
$PATH_VMATCH/mkvtree -db $REF_FA -v -pl -sti1 -bwt -dna -bck -suf -lcp -tis -ois -skp

## Using a length of the rep elements of insert size -20% and edit distance 5
cd ..
REP_LENGTH=$(( $INSERT_SIZE - ($INSERT_SIZE / 5) ))
OUT_NAME_BASE_VCF=${OUT_NAME_BASE}_LEN${REP_LENGTH}_ED5_vmatch
$PATH_VMATCH/vmatch -d -p -l $REP_LENGTH -e 1 -showdesc 30 mkvtree/$REF_IDX_NAME >${OUT_NAME_BASE_VCF}.ed1.txt
mv ${OUT_NAME_BASE_VCF}.ed1.txt ${OUT_NAME_BASE_VCF}.txt

## Convert vmatch .txt file to pseudo vcf file for lookup within vcf file
grep -v '#' ${OUT_NAME_BASE_VCF}.txt | awk '{print $2"\t"$3"\t"$1"\t.\t.\t.\t.\tCHROM="$6";LENGTH="$1";SVPOS="$7"\t."}' >200_1.txt
grep -v '#' ${OUT_NAME_BASE_VCF}.txt | awk '{print $6"\t"$7"\t"$5"\t.\t.\t.\t.\tCHROM="$2";LENGTH="$5";SVPOS="$3"\t."}' >200_2.txt
cat 200_1.txt 200_2.txt | sort -nk 2 >${OUT_NAME_BASE_VCF}.vcf
rm 200_1.txt 200_2.txt

## Compare vcf positions if they lie within repeat positions of vmatch.vcf file and tag them
file=../${BAM_NAME_BASE}_cov.vcf
FILE_BASE=$( basename $file)
FILENAME=${FILE_BASE%.vcf}
grep -e '^#' $file >../${FILENAME}_rep.vcf
perl $PATH_TO_VARIANTS_PERL/vcf_contains_pos2.pl $file ${OUT_NAME_BASE_VCF}.vcf >>../${FILENAME}_rep.vcf

## Find homopolymer stretches
cd ..
file=${FILENAME}_rep.vcf
FILE_BASE=$( basename $file )
FILENAME=${FILE_BASE%.vcf}
perl $PATH_TO_VARIANTS_PERL/filter_homopolymers.pl $file $REF_FA >${FILENAME}_hopo.vcf

###
###
### Search and tag snp accumulating regions (SARs)
###
###
SNPCOUNTMAX=4
SNPREGION=$(( $INSERT_SIZE * 2 ))
file=${FILENAME}_hopo.vcf
FILE_BASE=$( basename $file)
FILENAME=${FILE_BASE%.vcf}
perl $PATH_TO_VARIANTS_PERL/filter_multi_snps_2vcf_2.pl $file $SNPCOUNTMAX $SNPREGION >${FILENAME}_sar.vcf

###
###
### Count callers per position (CV2=2 callers support the variant at this position) / incl. unprecise structural variants (CSV2)
### add tags for BP to closeby variants
###
CPV=2
file=${FILENAME}_sar.vcf
FILE_BASE=$( basename $file)
FILENAME=${FILE_BASE%.vcf}
perl $PATH_TO_VARIANTS_PERL/filter_caller2pos.pl $file $CPV >${FILENAME}_cpv.vcf

###
###
### Apply MRA filter and tagging
###
###
file=${FILENAME}_cpv.vcf
FILE_BASE=$( basename $file)
FILENAME=${FILE_BASE%.vcf}
perl $PATH_TO_VARIANTS_PERL/filter_vcfs_2vcf.pl ${FILENAME}.vcf $MAA $MRA >${FILENAME}_tagged_${MRA}.vcf

### Filter tagged
TAG_FILTER="REP|HOP|CSV1|CV1|CNSV1"
file=${FILENAME}_tagged_${MRA}.vcf
bash $PATH_TO_VARIANTS_PERL/D03_filter_tags.sh $file "$TAG_FILTER" >${FILENAME}_filter_${MRA}.vcf



###
###
### Generate statistics
###
###
cat ${FILENAME}_filter_${MRA}.vcf | \
grep -e 'SNP' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>cov_snp_2.txt

## get Indel/SV coverage
cat ${FILENAME}_filter_${MRA}.vcf | \
grep -E 'DEL|INS|IND' | grep -E 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>cov_indelsmall_2.txt

## get SV coverage
cat ${FILENAME}_filter_${MRA}.vcf | \
grep -E 'DEL|INS|INV|DUP|ITX|CTX|LI|COMPLEX' | grep -Ev 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>cov_sv_2.txt

## get BP coverage
cat ${FILENAME}_filter_${MRA}.vcf | \
grep -E 'SVTYPE=BP' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>cov_bp_2.txt

echo "totcov <- read.table(\""${BAM_NAME_BASE}"_cov_total.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >cov_var.R
if [ $( cat cov_snp_2.txt | wc -l ) -gt 0 ]; then
  echo "snp <- read.table(\"cov_snp_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "snp <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat cov_indelsmall_2.txt | wc -l ) -gt 0 ]; then
  echo "indel <- read.table(\"cov_indelsmall_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "indel <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat cov_sv_2.txt | wc -l ) -gt 0 ]; then
  echo "sv <- read.table(\"cov_sv_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
else
  echo "sv <- matrix(c(0,0,0,0,0),ncol=5,byrow=TRUE)" >>cov_var.R
fi
if [ $( cat cov_bp_2.txt | wc -l ) -gt 0 ]; then
  echo "bp <- read.table(\"cov_bp_2.txt\", sep=\"\t\", header=FALSE, dec=\".\")" >>cov_var.R
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

## plot coverage
echo "plot(subtot[c(2,3)],col=\"lightblue\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.3);" >>cov_var.R
echo "if( val %in% mysum ){ points(sub1[c(2,4)],col=\"red\",main=val,xlab=\"position\",ylab=\"coverage\",cex=.6); mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% indelsum ){ points(subindel[c(2,4)],col=\"blue\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100));mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% svsum ){ points(subsv[c(2,4)],col=\"orange\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% bpsum ){ points(subbp[c(2,4)],col=\"black\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "par(xpd=TRUE)" >>cov_var.R
echo "legend(\"bottomleft\", inset=c(-0.18,-0.6), c(\"SNP\",\"InDel\",\"SV\",\"BP\"),lty=c(NA,NA,NA,NA),lwd=c(2.5,2.5,2.5,2.5),col=c(\"red\",\"blue\",\"orange\",\"black\"),pch=c(1,1,1,1),cex=.6)" >>cov_var.R
echo "par(xpd=FALSE)" >>cov_var.R

## plot frequency
echo "if( val %in% mysum ){ plot(sub1[c(2,5)],col=\"red\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100));mtext(fname,side=1,line=4,cex=.6) }" >>cov_var.R
echo "if( val %in% indelsum ){ points(subindel[c(2,5)],col=\"blue\",main=val,xlab=\"position\",ylab=\"frequency\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% svsum ){ points(subsv[c(2,5)],col=\"orange\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "if( val %in% bpsum ){ points(subbp[c(2,5)],col=\"black\",cex=.6,ylim=c(0,100)) }" >>cov_var.R
echo "}" >>cov_var.R

Rscript cov_var.R
mv Rplots.pdf "${BAM_NAME_BASE}"_coverages.pdf

## plot variant frequencies
less ${FILENAME}_filter_${MRA}.vcf | \
grep -e SNP | grep -Ev '^#|REP|CV1' | grep -Ev 'samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[3]"\t"pct[4] }' | sort -un -k2,2 \
>${BAM_NAME_BASE}_pct_hist.txt
sed -i '1ichrom\tpos\tcov\tfreq' ${BAM_NAME_BASE}_pct_hist.txt

## generate file for rscript
echo "snp <- read.table(\""${BAM_NAME_BASE}_pct_hist.txt"\", sep=\"\t\", header=TRUE, dec=\".\")" >pct_hist.R
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
mv Rplots.pdf ${BAM_NAME_BASE}_frequency.pdf
