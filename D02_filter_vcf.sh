#!/bin/bash

# module load vmatch

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# export variables for the use in perl
export PATH_CALLER_COLLECT=${PATH_CALLER_COLLECT}
export BAM_NAME_BASE=${BAM_NAME_BASE}
export REPEATS=${REPEATS}
export MIN_CPC=${MIN_CPC}
export MIN_CPV=${MIN_CPV}
export REF_FA=${REF_FA}


### runs collect_variants_varcap.pl to collect variants from the different callers into one vcf
PATH_TO_VARIANTS_PERL=$PATH_SCRIPTS/vcffilter
FILENAME_BASE=${FILENAME_BASE}_v1
echo $FILENAME_BASE

cd $PATH_PROJECTS_DATA/$PROJ_NAME
echo "VARCAP: collect variants"
mkdir -p vcfs
mkdir -p vcfs_raw
mkdir -p vcfs_raw/vcfs_temp
rm vcfs/*
rm vcfs_raw/*
rm vcfs_raw/vcfs_temp/*

### 1. Get raw vcfs/output files of all callers and join them into a master vcf
cd vcfs_raw
perl -I $PATH_TO_VARIANTS_PERL $PATH_TO_VARIANTS_PERL/collect_variants_varcap_2vcf_02.pl $PATH_VARCAP $PATH_CALLER $FILENAME_BASE >$PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/$BAM_NAME_BASE.info
# 1.1 unify chrom names,as different tools handle chrom names differently

# convert cortex chrom to ref chrom name
# cat "${BAM_NAME_BASE}_[0-9].vcf" | grep -v '^#' | cut -f1,8 | sort -u -k1,1 >chrom_table.txt
# cat chrom_table.txt | grep -e 'cortex' >chrom_table_cortex.txt
>${BAM_NAME_BASE}_${COUNTER}.vcf.mod
COUNTER=1
while  [ ${COUNTER} -le ${REPEATS} ]; do
  CHROMS=$( cat "${BAM_NAME_BASE}_${COUNTER}.vcf" | grep -Ev '^#|cortex' | cut -f1 | sort -u )
  cat "${BAM_NAME_BASE}_${COUNTER}.vcf" | while read -r line; do
    if [[ $line == *cortex* ]]; then
      CCHR=$( echo -e "$line" | cut -f1 )
      CREST=$( echo -e "$line" | cut -f2- )
      NCHR=$( echo -e "$CHROMS" | grep -e "$CCHR" )
      echo -e "$NCHR\t$CREST" >>${BAM_NAME_BASE}_${COUNTER}.vcf.mod
    else
      echo -e "$line" >>${BAM_NAME_BASE}_${COUNTER}.vcf.mod
    fi
  done
  # sort according to chrom and pos
  >HEAD.txt
  >BODY.txt
  cat ${BAM_NAME_BASE}_${COUNTER}.vcf.mod | while read -r line; do
    if [[ $line ==  \#* ]]; then
      echo -e "$line" >>HEAD.txt
    else
      echo -e "$line" >>BODY.txt
    fi
  done
  cat BODY.txt | sort -k1,1 -k2,2n >BODY.sort.txt
  cat HEAD.txt BODY.sort.txt >${BAM_NAME_BASE}_${COUNTER}.vcf.mod
  rm HEAD.txt BODY.sort.txt BODY.txt
  mv ${BAM_NAME_BASE}_${COUNTER}.vcf.mod ${BAM_NAME_BASE}_${COUNTER}.vcf
  let COUNTER=$COUNTER+1
done


# # ${BAM_NAME_BASE}_${COUNTER}.vcf
# COUNTER=1
# while  [ ${COUNTER} -le ${REPEATS} ]; do
  # iterate through replicated
#  while read line; do
#	if [[ $line == \#* ]]; then
#	  echo -e "$line"
#	else
#	  CHROM1=$( echo -e "$line" | cut -f1 )
#	  REST1=$( echo -e "$line" | cut -f2- )
#	  if [[ $CHROM1 == gi\|* ]]; then
#	    CHROM2=$( echo -e "$CHROM1" | cut -d '|' -f 4  )
#	    echo -e "$CHROM2\t$REST1"
#	  else
#	    echo -e "$CHROM1\t$REST1"
#	  fi
#	fi
  
#  done < "${BAM_NAME_BASE}_${COUNTER}.vcf" >${BAM_NAME_BASE}_${COUNTER}.vcf.upd
  # mv ${BAM_NAME_BASE}_${COUNTER}.vcf.upd ${BAM_NAME_BASE}_${COUNTER}.vcf
# let COUNTER=$COUNTER+1
# done


### 2. Calculate coverages and percentages for positions and add them to vcf
COUNTER=1
while  [ ${COUNTER} -le ${REPEATS} ]; do
  echo "VARCAP run $COUNTER: calculate average coverage per chromosome/contig"
  $PATH_SAMTOOLS/samtools depth $PATH_BWA_DATA/${BAM_NAME_BASE}_bwa_${COUNTER}.bam | awk 'NR%500==0' >${BAM_NAME_BASE}_${COUNTER}_cov_total.txt
  # CHROM=$( cat ${BAM_NAME_BASE}_${COUNTER}_cov_total.txt | cut -f1 | sort -u )
  cat "${BAM_NAME_BASE}_${COUNTER}_cov_total.txt" | awk '{ a[$1]+=$3; b[$1]=NR; next } END {for(i in a) { av_cov=(a[i]/b[i]); print i"\t"av_cov } }' >${BAM_NAME_BASE}_${COUNTER}_cov_av.txt

  echo "VARCAP run $COUNTER: update total coverage for variant positions"
  echo "${COUNTER}:get coverages for positions"
  grep -v '#' ${BAM_NAME_BASE}_*${COUNTER}.vcf | cut -f 1,2 >pos_${COUNTER}.bed
  $PATH_SAMTOOLS/samtools depth -b pos_${COUNTER}.bed $PATH_BWA_DATA/${BAM_NAME_BASE}_bwa_${COUNTER}.bam >${BAM_NAME_BASE}_${COUNTER}_cov_pos.txt
  echo "${BAM_NAME_BASE}_${COUNTER}: update file"
  perl -I $PATH_TO_VARIANTS_PERL $PATH_TO_VARIANTS_PERL/get_coverage_2vcf2.pl "${BAM_NAME_BASE}_${COUNTER}.vcf" "${BAM_NAME_BASE}_${COUNTER}_cov_pos.txt" "${BAM_NAME_BASE}_${COUNTER}_cov_av.txt" >${BAM_NAME_BASE}_bwa_${COUNTER}_cov.vcf
  cp ${BAM_NAME_BASE}_bwa_${COUNTER}_cov.vcf vcfs_temp/${BAM_NAME_BASE}_bwa_${COUNTER}_cov.vcf
  cp ${BAM_NAME_BASE}_bwa_${COUNTER}.vcf ${BAM_NAME_BASE}_bwa_${COUNTER}_cov.vcf ../vcfs
  let COUNTER=$COUNTER+1
done

### 3. Search for repetitive elements within the reference genome that are longer than insert size and tag homopolymers
echo "VARCAP: search for repetitive elements using vmatch"
mkdir -p vmatch/mkvtree
cd vmatch/mkvtree
REF_IDX_NAME=$( basename $REF_FA )
OUT_NAME_BASE=$( echo $REF_IDX_NAME | sed 's/\..*$//' )
mkvtree -db $REF_FA -v -pl -sti1 -bwt -dna -bck -suf -lcp -tis -ois -skp

# vmatch using a length of the rep elements of insert size -20% and edit distance 5
cd ..
REP_LENGTH=$(( $INSERT_SIZE - ($INSERT_SIZE / 5) ))
echo "REP_LENGTH="$REP_LENGTH
OUT_NAME_BASE_VCF=${OUT_NAME_BASE}_LEN${REP_LENGTH}_ED5_vmatch
vmatch -d -p -l $REP_LENGTH -e 1 -showdesc 30 mkvtree/$REF_IDX_NAME >${OUT_NAME_BASE_VCF}.ed1.txt
mv ${OUT_NAME_BASE_VCF}.ed1.txt ${OUT_NAME_BASE_VCF}.txt

# convert vmatch .txt file to pseudo vcf file for lookup within vcf file
grep -v '#' ${OUT_NAME_BASE_VCF}.txt | awk '{print $2"\t"$3"\t"$1"\t.\t.\t.\t.\tCHROM="$6";LENGTH="$1";SVPOS="$7"\t."}' >200_1.txt
grep -v '#' ${OUT_NAME_BASE_VCF}.txt | awk '{print $6"\t"$7"\t"$5"\t.\t.\t.\t.\tCHROM="$2";LENGTH="$5";SVPOS="$3"\t."}' >200_2.txt
cat 200_1.txt 200_2.txt | sort -nk 2 >${OUT_NAME_BASE_VCF}.vcf
rm 200_1.txt 200_2.txt

# compare vcf positions if they lie within repeat positions of vmatch.vcf file and tag them
echo "VARCAP: search and tag repetitive regions within vcf file"
COUNTER=1
for file in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf ); do
  # first get header of vcf file
  FILE_BASE=$( basename $file)
  FILENAME=${FILE_BASE%.vcf}
  grep -e '^#' $file >../${FILENAME}_rep.vcf
  # then run/add calls
  perl ${PATH_SCRIPTS}/vcffilter/vcf_contains_pos2.pl $file ${OUT_NAME_BASE_VCF}.vcf >>../${FILENAME}_rep.vcf
  rm $file
  cp ../${FILENAME}_rep.vcf ../vcfs_temp/
  let COUNTER=$COUNTER+1
done

# find homopolymer stretches
echo "VARCAP: search and tag homopolymer regions within vcf file"
cd $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw
COUNTER=1
for file2 in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf ); do
  # first get header of vcf file
  FILE_BASE=$( basename $file2 )
  FILENAME=${FILE_BASE%.vcf}
  perl $PATH_SCRIPTS/vcffilter/filter_homopolymers.pl $file2 $REF_FA >${FILENAME}_hopo.vcf
  rm $file2
  cp ${FILENAME}_hopo.vcf vcfs_temp/
done

### 4. Search and tag snp accumulating regions (SARs)
SNPCOUNTMAX=4
SNPREGION=$(( $INSERT_SIZE * 2 ))

cd $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw
echo "VARCAP: search and tag snp accumulating regions (SARs)"
COUNTER=1
for file in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf ); do
  FILE_BASE=$( basename $file)
  FILENAME=${FILE_BASE%.vcf}
  perl ${PATH_SCRIPTS}/vcffilter/filter_multi_snps_2vcf_2.pl $file $SNPCOUNTMAX $SNPREGION >${FILENAME}_sar.vcf
  rm $file
  cp ${FILENAME}_sar.vcf vcfs_temp/
  # echo "$COUNTER"
  let COUNTER=$COUNTER+1
done

# add tags for BP to closeby variants, count callers per position (CV2=2 callers support the variant at this position) / incl. unprecise structural variants (CSV2)
echo "Mark callers/position (CV[1..n]) incl. unprecise SV(CSV[1..n]) and BP accompanying variants (BPA)"

COUNTER=1
CPV=2 # caller per variant
for file in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf ); do
  FILE_BASE=$( basename $file)
  FILENAME=${FILE_BASE%.vcf}
  perl ${PATH_SCRIPTS}/vcffilter/filter_caller2pos.pl $file $CPV >${FILENAME}_cpv.vcf
  rm $file
  cp ${FILENAME}_cpv.vcf vcfs_temp/
  # echo "$COUNTER"
  let COUNTER=$COUNTER+1
done

### 5. Apply MRA filter
echo "VARCAP: filter vcf file"

if [ "$#" -gt 0 ];then
  MRA=(${@:1})
  echo "VARCAP: Use input MRA: ${MRA[*]}"
  # exit
else
  MRA=($( cat $PATH_PROJECTS_DATA/$PROJ_NAME/variant.config | grep -e '^MRA=' | cut -d'=' -f2 ))
  echo "VARCAP: Use default MRA: ${MRA[*]}"
fi
MAA=$( cat $PATH_PROJECTS_DATA/$PROJ_NAME/variant.config | grep -e '^MAA=' | cut -d'=' -f2 )
CUTOFF_AB=$MAA

COUNTER=1
# loop through files
for file in $( ls $PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/vcfs_temp/*.vcf ); do
  FILE_BASE=$( basename $file)
  FILENAME=${FILE_BASE%.vcf}

  # loop through cutoffs for read coverage
  for CUTOFF_PC in "${MRA[*]}"
  do
    echo -e "MAA:$CUTOFF_AB MRA:$CUTOFF_PC"
    perl $PATH_SCRIPTS/vcffilter/filter_vcfs_2vcf.pl ${FILENAME}.vcf $CUTOFF_AB $CUTOFF_PC >${FILENAME}_filter_${CUTOFF_PC}.vcf
    perl $PATH_SCRIPTS/vcffilter/filter_vcfs_2vcf.pl ${FILENAME}.vcf $CUTOFF_AB 0 >${FILENAME}_filter_none.vcf
    # filter acccording to tags
    TAG_FILTER="REP|HOP|CSV1|CV1|CNSV1"
    bash $PATH_SCRIPTS/../D03_filter_tags.sh ${FILENAME}_filter_${CUTOFF_PC}.vcf "$TAG_FILTER" >${FILENAME}_filter_${CUTOFF_PC}tags.vcf
    
    cp ${FILENAME}_filter_${CUTOFF_PC}tags.vcf vcfs_temp/${FILENAME}_filter_${CUTOFF_PC}tags.vcf
    cp ${FILENAME}_filter_${CUTOFF_PC}.vcf ../vcfs
    cp ${FILENAME}_filter_${CUTOFF_PC}tags.vcf ../vcfs
  done
  rm $file
  let COUNTER=$COUNTER+1
done

# 5.1 unify chrom names,as different tools handle chrom names differently
# ${BAM_NAME_BASE}_${COUNTER}.vcf
COUNTER=1
while  [ ${COUNTER} -le ${REPEATS} ]; do
  # iterate through replicated
  while read line; do
        if [[ $line == \#* ]]; then
          echo -e "$line"
        else
          CHROM1=$( echo -e "$line" | cut -f1 )
          REST1=$( echo -e "$line" | cut -f2- )
          if [[ $CHROM1 == gi\|* ]]; then
            CHROM2=$( echo -e "$CHROM1" | cut -d '|' -f 4  )
            echo -e "$CHROM2\t$REST1"
          else
            echo -e "$CHROM1\t$REST1"
          fi
        fi

  done < "${FILENAME}_filter_${CUTOFF_PC}.vcf" >${FILENAME}_filter_${CUTOFF_PC}_upd.vcf
  mv ${FILENAME}_filter_${CUTOFF_PC}_upd.vcf ${FILENAME}_filter_${CUTOFF_PC}.vcf
  cp ${FILENAME}_filter_${CUTOFF_PC}.vcf ../vcfs
let COUNTER=$COUNTER+1
done

### 6. Generate statistics
### 6.1 print coverage
# get genome coverage
# coverages are stored in vcfs_raw/*1_cov_total.txt

# get coverage avg and stdev of variants
# get SNP coverage
cd $PATH_PROJECTS_DATA/$PROJ_NAME
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*tags.vcf | \
grep -e 'SNP' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_snp_2.txt
# get Indel/SV coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*tags.vcf | \
grep -E 'DEL|INS|IND' | grep -E 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | grep -Ev '^#|samtools|gatk|cortex' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_indelsmall_2.txt
# get SV coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*tags.vcf | \
grep -E 'DEL|INS|INV|DUP|ITX|CTX|LI|COMPLEX' | grep -Ev 'SVLEN=[-]{0,1}[0-9]{1}\W' | grep -Ev 'REP|CV1' | awk '{ split($10, pct, ":"); print $1"\t"$2"\t"pct[1]"\t"pct[3]"\t"pct[4] }' \
>vcfs_raw/cov_sv_2.txt
# get BP coverage
cat vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*tags.vcf | \
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
less vcfs/${BAM_NAME_BASE}_bwa_1_*_filter_*tags.vcf | \
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
