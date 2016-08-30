#!/bin/bash

BAM=$1
BAM_BASE=$( basename $BAM | sed 's/.bam$//' )

PATH_CORTEX_DATA=
REF=
REF_NAME=$(basename $REF | sed 's/\.f.*a$//')
PATH_PICARD=/opt/apps/picard_tools/2.0.1
PATH_VCFTOOLS=/opt/apps/vcftools/vcftools_0.1.9
PATH_STAMPY=/opt/apps/stampy/stampy-1.0.21
PATH_CORTEX_DIR=/opt/apps/cortex/1.0.5.14
PATH_CORTEX=/opt/apps/cortex/1.0.5.14/bin
PATH_RUN_CALLS=/opt/apps/cortex/1.0.5.14/scripts/calling

###
READS1=${BAM_BASE}_1.fastq
READ_NAME=$(basename $READS1 | sed 's/..\.f.*q$//')
READS2=${BAM_BASE}_2.fastq

# modify libs
PERL5LIB=$PATH_CORTEX_DIR/scripts/analyse_variants/bioinf-perl/lib:$PATH_CORTEX_DIR/scripts/calling:$PERL5LIB
PATH=$PATH_CORTEX_DIR/scripts/analyse_variants/needleman_wunsch:$PATH

# 0. prepare files
echo "Create cortex dir: $PATH_CORTEX_DATA"
mkdir -p $PATH_CORTEX_DATA

mkdir -p $PATH_CORTEX_DATA/
mkdir -p $PATH_CORTEX_DATA/$BAM_BASE
cd $PATH_CORTEX_DATA
## replace/add read group header
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PATH_PICARD/picard.jar AddOrReplaceReadGroups I=$BAM O=${BAM_BASE}.rgroup.bam LB=D PL=illumina PU=S SM=B VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false
##convert bam file to fastq reads
java -jar -Xmx3g $PATH_PICARD/picard.jar SamToFastq INPUT=${BAM_BASE}.rgroup.bam FASTQ=${READS1} SECOND_END_FASTQ=${READS2} INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT

# 1. prepare cortex index, if not exists
if [[ ! -e $PATH_CORTEX_DATA/ref/ref.k31.ctx ]]; then
  mkdir -p $PATH_CORTEX_DATA
  cd $PATH_CORTEX_DATA
  rm $PATH_CORTEX_DATA/*
  rm $PATH_CORTEX_DATA/ref/*
  ### prepare
  # this script sets up the preriquisites for running cortex with the given reference
  
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

fi

# 2. modify input files file
cd $PATH_CORTEX_DATA
sed -i 's#^.*$#'$READS1'#' pe_filelist1
sed -i 's#^.*$#'$READS2'#' pe_filelist2

# 3. run run_calls.pl
perl  $PATH_RUN_CALLS/run_calls.pl --first_kmer 31 --last_kmer 51 --kmer_step 20 --fastaq_index ${REF_NAME}_index --auto_cleaning yes --bc yes --pd no --outdir $PATH_CORTEX_DATA/$BAM_BASE --outvcf $READ_NAME \
--ploidy 1 --stampy_hash $REF_NAME --stampy_bin $PATH_STAMPY/stampy.py --list_ref_fasta file_listing_fasta --refbindir ref/ \
--genome_size 2414465 --max_read_len 10000 --qthresh 5 --mem_height 21 --mem_width 100 --vcftools_dir $PATH_VCFTOOLS \
--do_union yes --ref CoordinatesAndInCalling --workflow independent --logfile log.txt,f


rm ${BAM_BASE}.rgroup.bam
