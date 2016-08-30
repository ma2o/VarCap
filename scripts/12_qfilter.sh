#!/bin/bash
#$-q all.q@cube[ab]*
#$-N filterFastq
#$-l vf=10G
#$-o log
#$-e log

INFILE=$1
OUTDIR_QC=fastqc
OUTDIR_PRINSEQ=prinseq_data
TRIM_R=$2
OUTFILE=$OUTDIR_PRINSEQ/$(basename ${INFILE} | sed 's/\.f.*$//')
OUTFILE_QC=$(basename ${INFILE} | sed 's/\.f.*$//')

mkdir $OUTDIR_QC
mkdir $OUTDIR_PRINSEQ

#conservative filter
prinseq-lite -fastq $INFILE -trim_right $TRIM_R -out_good $OUTFILE"_trr" -out_bad $OUTFILE"_bad_trr" >$OUTFILE"_info_trr".txt
prinseq-lite -fastq $OUTFILE"_trr".fastq -trim_qual_right 20 -out_good $OUTFILE"_trr_trq20" -out_bad $OUTFILE"_bad_trr_trq20" >$OUTFILE"_info_trr_trq20".txt
prinseq-lite -fastq $OUTFILE"_trr_trq20".fastq -min_qual_mean 30 -out_good $OUTFILE_QC"_filter" -out_bad $OUTFILE"_filter_bad" >$OUTFILE"_filter_info".txt
fastqc -o $OUTDIR_QC/ $OUTFILE_QC"_filter".fastq
