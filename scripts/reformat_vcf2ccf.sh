#!/bin/bash
#$-q all.q@cube[ab]*

PROJECT_ID=$1
EXPERIMENT_ID=$2
SEQUENCING_ID=$3
EXPERIMENT=$4
REPLICATE=$5
SELECTION=$6
SELECTION_VALUE=$7
TIMEPOINT=$8
FILE_NAME=$9
EXP_ITERATION=${10}
FILE1=${11}
PATH_TO_SCRIPT=${12}

perl ${PATH_TO_SCRIPT}/convert_vcf2mysql2.pl $FILE1 | awk -v exp_id=$EXPERIMENT_ID -v iter=$EXP_ITERATION -v proj_id=$PROJECT_ID -v seq_id=$SEQUENCING_ID -v file_name=$FILE_NAME -v experiment=$EXPERIMENT -v replicate=$REPLICATE -v selection=$SELECTION -v selection_value=$SELECTION_VALUE -v timepoint=$TIMEPOINT '{print proj_id"\t"exp_id"\t"seq_id"\t"experiment"\t"replicate"\t"selection"\t"selection_value"\t"timepoint"\t"file_name"\t"iter"\t"$0}'


