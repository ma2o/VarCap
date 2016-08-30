#!/bin/bash
#$-q all.q@cube[ab]*

BAM=$1
BAM_NAME=$(basename $BAM | sed 's/\..am$//')
OUTDIR=$2
REF_FA=$3
REF_FA_BASE=$( basename $REF_FA )
PATH_DELLY=$4
PATH_DELLY_DATA=$5
INSERT_SIZE=$6
REF_FA_MOD=$7
INSERT_SIZE_CUTOFF=$( echo $(( $INSERT_SIZE + $INSERT_SIZE / 3 )) )
REF=$REF_FA_MOD/$REF_FA_BASE
#
cd $PATH_DELLY_DATA/$OUTDIR
# from -q 20 to -q 0
$PATH_DELLY/delly -o $BAM_NAME.del.txt -b $BAM_NAME.br.txt -q 0 -s $INSERT_SIZE_CUTOFF -g $REF -p $BAM

# $PATH_DELLY/duppy -o $BAM_NAME.dup.txt -b $BAM_NAME.dup.br.txt -q 20 -g $REF -p $BAM
$PATH_DELLY/duppy -o $BAM_NAME.dup.txt -b $BAM_NAME.dup.br.txt -q 0 -g $REF -m 13 -p $BAM
# from -q 10 to 0
$PATH_DELLY/invy -o $BAM_NAME.inv.txt -r $BAM_NAME.inv_merge.txt -b $BAM_NAME.inv.br.txt -k $BAM_NAME.inv.br_merge.txt -q 0 -g $REF -p $BAM
# from -q 15 to 0
$PATH_DELLY/jumpy -o $BAM_NAME.jmp.txt -r $BAM_NAME.jmp_merge.txt -b $BAM_NAME.jmp.br.txt -k $BAM_NAME.jmp.br_merge.txt -q 0 -g $REF -p $BAM


