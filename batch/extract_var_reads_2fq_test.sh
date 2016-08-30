#!/bin/bash

BAM=$1
VCF=$2
REF=$3
ALLIN=( "$@" )
PCT=( ${ALLIN[@]:3:${#ALLIN[@]}} )
# input format: bam vcf ref lo1 up1 lo2 up2

echo "$BAM $VCF ${PCT[*]}"
# create output names/dir
BAM_NAME=$( basename $BAM | sed 's/\.bam$//' )

# loop over frequencies and create fq.gz files
PCT_TEMP=( ${PCT[@]} )
while [ ${#PCT_TEMP[@]} -gt "1" ]; do
  PCT_EX=(${PCT_TEMP[@]:0:2})
  PCT_RE=(${PCT_TEMP[@]:2})
  PCT_TEMP=(${PCT_RE[@]})
  echo "pctex:${PCT_EX[@]}"
  # echo "pctte:${PCT_TEMP[@]}"
  # calculate average of pct value, include it into name
  PCT_AV=$(( ( ${PCT_EX[0]} + ${PCT_EX[1]} ) / 2 ))
  # echo "pct_av $PCT_AV"
  # generate name
  OUTNAME=$( echo "${BAM_NAME}_${PCT_AV}" )
  echo "$OUTNAME"
  # generate frequency filtered vcf files 
  cat $VCF | grep -Ev '^#|cortex' | grep -e 'SNP' | grep -v 'CV1' | awk -v low=${PCT_EX[0]} -v high=${PCT_EX[1]} '{ split($10,abs,":"); freq=abs[4]; if (low <= freq && freq <= high) {print} }' >${OUTNAME}.vcf
  # extract reads based on vcf file positions and cigar strings of bam files
    >$OUTNAME.sam
	echo "Extract reads"
	cat ${OUTNAME}.vcf | while read -r line; do
		if [[ $line != '#'* ]]; then
			CHROM=$( echo $line | cut -d" " -f1 ) 
			POS=$( echo $line | cut -d" " -f2 )
			REF=$( echo $line | cut -d" " -f4 )
			ALT=$( echo $line | cut -d" " -f5 )
			# search pattern for cigar string
			PATTERN="[0-9]*[$ALT][0-9]*"
			# echo $PATTERN
			samtools view $BAM $CHROM:$POS-$POS | awk -v alt=$ALT -v pattern=$PATTERN -v pos=$POS '{ if ( $13 ~ /:[0-9]*[ACGT][0-9]*/ ) if ( $13 ~ pattern ) print $0 }' >>${OUTNAME}.sam
		fi
		
	done
    ### convert sam to fastq
    echo "Convert to and reheader bam file."
    samtools view -Sb -T $REF ${OUTNAME}.sam >${OUTNAME}.bam
    samtools fastq -t -1 ${OUTNAME}_1.fastq -2 ${OUTNAME}_2.fastq ${OUTNAME}.bam
    # remove temp files
    rm ${OUTNAME}.sam
    rm ${OUTNAME}.bam
done
