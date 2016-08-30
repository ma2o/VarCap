#!/bin/bash

# supply a list of sample names, which are used for further remap analysis
# correct format should be samplename with lower and upper boundaries for variant extraction
# samplename map_cov tot_snps var_snps lo1 hi1 lo2 hi2 ...

INFILE=$1
CURR_WD=$( pwd | sed 's/\/$//')

cat $INFILE | grep -v '^NAME' | while read -r line; do
# for file in $( cat $INFILE | cut -f1 | grep -v '^NAME' ); do 
  # FN=$( echo $file | sed 's/\/$//' )
  # sample name and percentages
  FN=$( echo -e "$line" | cut -f1 | sed 's/\/$//' )
  PCT=( $( echo -e "$line" | cut -f5- | tr '\t' ' ' ) )
  
  # create projects
  cd $CURR_WD/$FN
  # link ref
  REFNAME=$( cat variant.config | grep -e '^REF_FA_ORIGINAL' | cut -d'=' -f2 | sed 's#^.*\/##' )
  REFPATH=$CURR_WD/reference
  REF="$REFPATH/$REFNAME"
  echo "sample:$FN $REF"
  # create remap dir
  mkdir -p $CURR_WD/remap
  # create project directories +add percentage to samplename
  PCT_TEMP=( ${PCT[@]} )
  while [ ${#PCT_TEMP[@]} -gt "1" ]; do
	  PCT_EX=(${PCT_TEMP[@]:0:2})
	  PCT_RE=(${PCT_TEMP[@]:2})
	  PCT_TEMP=(${PCT_RE[@]})
	  PCT_AV=$(( ( ${PCT_EX[0]} + ${PCT_EX[1]} ) / 2 ))
	  OUTDIR=$( echo "${FN}_${PCT_AV}" )
	  echo "outdir:$OUTDIR"
	  mkdir $CURR_WD/remap/$OUTDIR
	  # create symbolic links to reference and bwa_index
	  ln -s $CURR_WD/reference reference
	  ln -s $CURR_WD/bwa_index bwa_index
	  # prepare dir and run refsearch
	  cp $CURR_WD/$FN/* $CURR_WD/remap/$OUTDIR/
	  mkdir -p $CURR_WD/remap/$OUTDIR/reference
	  mkdir -p $CURR_WD/remap/$OUTDIR/filter
	  mkdir -p $CURR_WD/remap/$OUTDIR/mapper/bwa
	  mkdir -p $CURR_WD/remap/$OUTDIR/logs
	  # change paths in variant.config
	  VPROJECT=$( cat $CURR_WD/remap/$OUTDIR/variant.config | grep -e '^PATH_PROJECTS_DATA' | cut -d'=' -f2 )
	  # PATH_PROJECTS_DATA
	  sed -i 's#^PATH_PROJECTS_DATA.*#PATH_PROJECTS_DATA='"$VPROJECT"'\/remap#' $CURR_WD/remap/$OUTDIR/variant.config
	  # PROJ_NAME
	  sed -i 's#^PROJ_NAME.*#PROJ_NAME='"$OUTDIR"'#' $CURR_WD/remap/$OUTDIR/variant.config
	  # BAM_NAME_BASE
	  sed -i 's#^BAM_NAME_BASE.*#BAM_NAME_BASE='"$OUTDIR"'#' $CURR_WD/remap/$OUTDIR/variant.config
	  # SUBSAMPLE_SIZE_ALT
	  sed -i 's#^SUBSAMPLE_SIZE_ALT.*#SUBSAMPLE_SIZE_ALT=-1#' $CURR_WD/remap/$OUTDIR/variant.config
	  # move to filter dir and extract reads
	  
	  cd $CURR_WD/remap/$OUTDIR/filter
	  echo "subscript extract reads."
	  bash $CURR_WD/extract_var_reads_2fq_test.sh $CURR_WD/$FN/mapper/bwa/${FN}_bwa_1.bam $CURR_WD/$FN/vcfs_raw/*filter_none.vcf $REF ${PCT_EX[0]} ${PCT_EX[1]}
	  mv *_1.fastq* ${FN}_${PCT_AV}_alt_1.fastq
          mv *_2.fastq* ${FN}_${PCT_AV}_alt_2.fastq
	  cd $CURR_WD/$FN
  done
  
done
