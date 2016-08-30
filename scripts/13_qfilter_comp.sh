#!/bin/bash

# SLURM
#SBATCH --job-name=VCqfilter
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=log/qfilter-%A_%a.out
#SBATCH --error=log/qfilter-%A_%a.err

INFILE=$1
OUTDIR_QC=fastqc
OUTDIR_PRINSEQ=prinseq_data
TRIM_R=$2
TRIM_L=$3
QUAL_WINDOW=$4
MIN_LENGTH=$5
PATH_PRINSEQ=$6
PATH_FASTQC=$7
OUTSUFF=$8
if [[ ! -z $OUTSUFF ]]; then
  OUTSUFF=alt_1
fi

# usage: bash 13_qfilter_comp.sh ${PROJ_NAME}_1.fastq $READS1_TRIM $READS1_TRIMF $QUAL_WINDOW_1 $MIN_LENGTH_1

OUTFILE_NAME=$(basename ${INFILE} | sed 's/\.f.*$//')

mkdir -p $OUTDIR_QC
mkdir -p $OUTDIR_PRINSEQ

## hard cutoff front and rear, quality window and overall length and quality
## check if variables for each filter steps are set, else skip this filtering step
# quality window from right
# if [ "$QUAL_WINDOW" -gt "0" ];
#   then
#     echo "trim qual window: $QUAL_WINDOW"
#     prinseq-lite -fastq $OUTFILE"_mod".fastq -trim_qual_window $QUAL_WINDOW -trim_qual_type min -trim_qual_right 20 -out_good $OUTFILE"_window" -out_bad $OUTFILE"_bad_win10" >$OUTFILE"_info_win10".txt
#     mv $OUTFILE"_window".fastq $OUTFILE"_mod".fastq
# fi
# minimum length filter
# if [ "$MIN_LENGTH" -gt "0" ];
#   then
#     echo "filter minimum length: $MIN_LENGTH"
#     prinseq-lite -fastq $OUTFILE"_mod".fastq -min_len $MIN_LENGTH -out_good $OUTFILE"_minlength" -out_bad $OUTFILE"_bad_minlength" >$OUTFILE"_info_minlength".txt
#     mv $OUTFILE"_minlength".fastq $OUTFILE"_mod".fastq
# fi

# concatenated prinseq filtering commands
echo "filtering reads trim_right:$TRIM_R trim_left:$TRIM_L window_size:$QUAL_WINDOW min_window_qual:20 min_length:$MIN_LENGTH min_qual_mean:30"
$PATH_PRINSEQ/prinseq-lite.pl -fastq $INFILE -trim_right $TRIM_R -trim_left $TRIM_L -trim_qual_window $QUAL_WINDOW -trim_qual_type min -trim_qual_right 20 -min_len $MIN_LENGTH -min_qual_mean 30 -out_good ${OUTFILE_NAME}"_"${OUTSUFF} -out_bad $OUTDIR_PRINSEQ/$OUTFILE_NAME"_bad" >$OUTDIR_PRINSEQ/$OUTFILE_NAME"_info_filter".txt

# run fastqc quality control
$PATH_FASTQC/fastqc -o $OUTDIR_QC/ $OUTFILE_NAME"_"${OUTSUFF}.fastq

