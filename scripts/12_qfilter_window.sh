#!/bin/bash
#$-q all.q@cube[ab]*
#$-N filterFastq
#$-l vf=10G
#$-o log
#$-e log

INFILE=$1
QUAL_WINDOW=$2
if [ ! $QUAL_WINDOW ];
  then
  QUAL_WINDOW="10"
fi
MIN_LENGTH=$3
if [ ! $MIN_LENGTH ];
  then
  MIN_LENGTH="35"
fi
OUTDIR_QC=fastqc
OUTDIR_PRINSEQ=prinseq_data
OUTFILE=$OUTDIR_PRINSEQ/$(basename ${INFILE} | sed 's/\.f.*$//')
OUTFILE_QC=$(basename ${INFILE} | sed 's/\.f.*$//' )

mkdir $OUTDIR_QC
mkdir $OUTDIR_PRINSEQ

#progressive filter
prinseq-lite -fastq $INFILE -trim_qual_window $QUAL_WINDOW -trim_qual_type min -trim_qual_right 20 -out_good $OUTFILE"_win10" -out_bad $OUTFILE"_bad_win10" >$OUTFILE"_info_win10".txt
prinseq-lite -fastq $OUTFILE"_win10".fastq  -min_len $MIN_LENGTH -out_good $OUTFILE"_win10_fi30" -out_bad $OUTFILE"_bad__win10_fi30" >$OUTFILE"_info_win10_fi30".txt
prinseq-lite -fastq $OUTFILE"_win10_fi30".fastq -min_qual_mean 30 -out_good $OUTFILE_QC"_filter" -out_bad $OUTFILE"_filter_bad" >$OUTFILE"_filter_info".txt
fastqc -o $OUTDIR_QC/ $OUTFILE_QC"_filter".fastq
