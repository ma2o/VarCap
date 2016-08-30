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
TRIM_L=$3
QUAL_WINDOW=$4
MIN_LENGTH=$5

# usage: bash 12_qfilter_front_rear_cut_window.sh ${PROJ_NAME}_1.fastq $READS1_TRIM $READS1_TRIMF $QUAL_WINDOW_1 $MIN_LENGTH_1

OUTFILE=$OUTDIR_PRINSEQ/$(basename ${INFILE} | sed 's/\.f.*$//')
OUTFILE_QC=$(basename ${INFILE} | sed 's/\.f.*$//')

mkdir $OUTDIR_QC
mkdir $OUTDIR_PRINSEQ

## hard cutoff front and rear, quality window and overall length and quality
## check if variables for each filter steps are set, else skip this filtering step
# trim right
if [ "$TRIM_R" -gt "0" ];
  then
    echo "trim right: $TRIM_R"
    prinseq-lite -fastq $INFILE -trim_right $TRIM_R -out_good $OUTFILE"_trimr" -out_bad $OUTFILE"_bad_trr" >$OUTFILE"_info_trr".txt
    mv $OUTFILE"_trimr".fastq $OUTFILE"_mod".fastq
  else
    mv $INFILE $OUTFILE"_mod".fastq
fi
# trim left
if [ "$TRIM_L" -gt "0" ];
  then
    echo "trim left: $TRIM_L"
    prinseq-lite -fastq $OUTFILE"_mod".fastq -trim_left $TRIM_L -out_good $OUTFILE"_triml" -out_bad $OUTFILE"_bad_trl" >$OUTFILE"_info_trl".txt
    mv $OUTFILE"_triml".fastq $OUTFILE"_mod".fastq
fi
# quality window from right
if [ "$QUAL_WINDOW" -gt "0" ];
  then
    echo "trim qual window: $QUAL_WINDOW"
    prinseq-lite -fastq $OUTFILE"_mod".fastq -trim_qual_window $QUAL_WINDOW -trim_qual_type min -trim_qual_right 20 -out_good $OUTFILE"_window" -out_bad $OUTFILE"_bad_win10" >$OUTFILE"_info_win10".txt
    mv $OUTFILE"_window".fastq $OUTFILE"_mod".fastq
fi
# minimum length filter
if [ "$MIN_LENGTH" -gt "0" ];
  then
    echo "filter minimum length: $MIN_LENGTH"
    prinseq-lite -fastq $OUTFILE"_mod".fastq -min_len $MIN_LENGTH -out_good $OUTFILE"_minlength" -out_bad $OUTFILE"_bad_minlength" >$OUTFILE"_info_minlength".txt
    mv $OUTFILE"_minlength".fastq $OUTFILE"_mod".fastq
fi
# minimum quality filter
echo "filter minimum quality mean: 30"
prinseq-lite -fastq $OUTFILE"_mod".fastq -min_qual_mean 30 -out_good $OUTFILE_QC"_filter" -out_bad $OUTFILE"_filter_bad" >$OUTFILE"_filter_info".txt

fastqc -o $OUTDIR_QC/ $OUTFILE_QC"_filter".fastq

