#!/bin/bash
#$-q all.q@cube[ab]*

. /etc/profile

READS1=$1
READ_NAME=$(basename $READS1 | sed 's/..\.f.*q$//')
READS2=$2
OUTDIR=$3
PATH_CORTEX_DATA=$4

cd $PATH_CORTEX_DATA
sed -i 's#^.*$#'$READS1'#' pe_filelist1
sed -i 's#^.*$#'$READS2'#' pe_filelist2
REF=$5
REF_NAME=$(basename $REF | sed 's/\.f.*a$//')
PATH_VCFTOOLS=$6
PATH_STAMPY=$7
PATH_CORTEX_DIR=$8
PATH_CORTEX=$9
PATH_RUN_CALLS=${10}

#add paths to PERL5LIB
PERL5LIB=$PATH_CORTEX_DIR/scripts/analyse_variants/bioinf-perl/lib:$PATH_CORTEX_DIR/scripts/calling:$PERL5LIB
PATH=$PATH_CORTEX_DIR/scripts/analyse_variants/needleman_wunsch:$PATH

cd $PATH_CORTEX_DATA
# 0: compile binaries in CORTEX_release_v1.0.5.21/bin
# make MAXK=63 NUM_COLS=2 cortex_var

#1: build reference genome binaries
# mkdir ref
# $PATH_CORTEX/cortex_var_31_c1 --kmer_size 31 --mem_height 17 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref.k31.ctx --sample_id REF
# $PATH_CORTEX/cortex_var_63_c1 --kmer_size 61 --mem_height 17 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref.k61.ctx --sample_id REF

# 2: build stampy hash of the genome
# $PATH_STAMPY/stampy.py -G $REF_NAME $REF
# $PATH_STAMPY/stampy.py -g $REF_NAME -H $REF_NAME

# build binary
# $PATH_CORTEX/cortex_var_31_c1 --pe_list pe_filelist1,pe_filelist2 --mem_height 17 --mem_width 100 --dump_binary test_bin.ctx --sample_id ZAM

# run cortex
# $PATH_CORTEX/cortex_var_63_c172

#3: run run_calls.pl
perl  $PATH_RUN_CALLS/run_calls.pl --first_kmer 31 --last_kmer 61 --kmer_step 30 --fastaq_index ${REF_NAME}_index --auto_cleaning yes --bc yes --pd no --outdir $OUTDIR --outvcf $READ_NAME \
--ploidy 1 --stampy_hash $REF_NAME --stampy_bin $PATH_STAMPY/stampy.py --list_ref_fasta file_listing_fasta --refbindir ref/ \
--genome_size 2414465 --max_read_len 10000 --qthresh 5 --mem_height 19 --mem_width 100 --vcftools_dir $PATH_VCFTOOLS \
--do_union yes --ref CoordinatesAndInCalling --workflow independent --logfile log.txt,f
#unused: --apply_pop_classifier --dups
