#!/bin/bash
#$-q all.q@cube[ab]*

. /etc/profile

# this script sets up the preriquisites for running cortex with the given reference
PATH_CORTEX_DATA=$1
REF=$2
REF_NAME=$(basename $REF | sed 's/\.f.*a$//')
PATH_VCFTOOLS=$3
PATH_STAMPY=$4
PATH_CORTEX_DIR=$5
PATH_CORTEX=$6
PATH_RUN_CALLS=$7

#add paths to PERL5LIB
PERL5LIB=$PATH_CORTEX_DIR/scripts/analyse_variants/bioinf-perl/lib:$PATH_CORTEX_DIR/scripts/calling:$PERL5LIB
PATH=$PATH_CORTEX_DIR/scripts/analyse_variants/needleman_wunsch:$PATH


cd $PATH_CORTEX_DATA

# NOTE: binaries have been precompiled, see below:
#0: compile binaries in CORTEX_release_v1.0.5.14/bin
#make MAXK=63 NUM_COLS=2 cortex_var

# create dummy starting files
echo "dummy_reads1" >pe_filelist1
echo "dummy_reads2" >pe_filelist2
echo -e  ${REF_NAME}"\t.\tpe_filelist1\tpe_filelist2" >${REF_NAME}_index
echo $REF >file_listing_fasta

#1: build reference genome binaries was ( --mem_height 17 )
mkdir ref
$PATH_CORTEX/cortex_var_31_c1 --kmer_size 31 --mem_height 21 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref/ref.k31.ctx --sample_id REF
$PATH_CORTEX/cortex_var_63_c1 --kmer_size 51 --mem_height 21 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref/ref.k61.ctx --sample_id REF

#2: build stampy hash of the reference genome
$PATH_STAMPY/stampy.py -G $REF_NAME $REF
$PATH_STAMPY/stampy.py -g $REF_NAME -H $REF_NAME

