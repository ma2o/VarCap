#!/bin/bash

# SLURM
#SBATCH --job-name=VCrefsearch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=reference/refsearch-%A_%a.out
#SBATCH --error=reference/refsearch-%A_%a.err

# read config file(variant.config) within the same directory
# CURRENT_DIR=$2
CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# path to references directory ( references in fasta format )
REFPATH=$1
# if defined, this references will be used for mapping and statistics, but excluded from reference determination
EXREFPATH=$2

# exclude those references from mapping
# EXCLUDE="human|mouse"
EXCLUDE="alien"

QUERY=$PATH_ALT_READS1
OUTDIR=$PATH_PROJECTS_DATA/$PROJ_NAME/reference
SUBSAMPLE_SIZE_TEST=100000
# check read depth of query reads, if lower than sub size adjust accordingly
QU_READS=$( zcat $QUERY | awk 'END {reads=NR/4; print reads }' )
if [[ $SUBSAMPLE_SIZE_TEST -gt $QU_READS ]]; then
  SUBSAMPLE_SIZE_TEST=$QU_READS
fi
OUTNAME=$( basename $QUERY | sed 's/\.fastq.*$//' )

CURR_DIR=$( pwd | sed 's/\/$//' )

mkdir -p $OUTDIR

# generate blacklist from references in EXREFPATH, empty if none specified
BLACKLIST=$( ls $EXREFPATH | tr "\t" "\n" | grep -E 'fasta$|fna$|fa$' | sed 's/\.f.*a$//' )
echo "$BLACKLIST" >$OUTDIR/blacklist.txt


echo "Inputs..."
echo "Ref:$REFPATH"
echo "Quy:$QUERY"
>$OUTDIR/$OUTNAME.info.txt

# generate subsample of query
if [ -f "$OUTDIR/$OUTNAME.fasta" ]; then
  echo "QUERY: $OUTNAME.fasta reference found, using it."
else
  echo "Subsampling..."
  python $PATH_SCRIPTS/filter/subsample.py -n $SUBSAMPLE_SIZE_TEST -in $QUERY -out $OUTDIR/$OUTNAME.fastq
  gzip $QUERY $QUERY.gz
  # fastq to fasta
  echo "Fastq to fasta..."
  python $PATH_SCRIPTS/filter/fastq2fasta.py < $OUTDIR/$OUTNAME.fastq >$OUTDIR/$OUTNAME.fasta
fi


map_bwa(){
  local REFPATH=$1
  local REF_IDX_NAME=$2
  mkdir -p $PATH_PROJECTS_DATA/bwa_index_ref
  # mkdir -p $PATH_PROJECTS_DATA/bwa_index/bwa_index_${REF_IDX_NAME}
  echo "searchRef_INPUT:$REFPATH $REF_IDX_NAME"
  echo "searchRef_IDXDIR:$PATH_PROJECTS_DATA/bwa_index/bwa_index_${REF_IDX_NAME}"
  # check if ref exists, else create it
  if [ -d "$PATH_PROJECTS_DATA/bwa_index/bwa_index_${REF_IDX_NAME}" ]; then
    echo "searchRef_REF: $REF_IDX_NAME index exists, using it."
  else
    
    mkdir -p $PATH_PROJECTS_DATA/bwa_index/bwa_index_${REF_IDX_NAME}
    $PATH_BWA_075/bwa index -p $PATH_PROJECTS_DATA/bwa_index/bwa_index_$REF_IDX_NAME/$REF_IDX_NAME $REFPATH
  fi
  
  # run bwa: add -x intractg for combined -B9 -O16 -L5  (intra-species contigs to ref) settings of bwa mem
  $PATH_BWA_075/bwa mem -x intractg -M -t 2 $PATH_PROJECTS_DATA/bwa_index/bwa_index_$REF_IDX_NAME/$REF_IDX_NAME $OUTDIR/$OUTNAME.fasta >$OUTDIR/$OUTNAME.sam
  # sam to bam
  samtools="$PATH_SAMTOOLS/samtools"
  $samtools view -bS $OUTDIR/${OUTNAME}.sam >$OUTDIR/$OUTNAME.bam
  $samtools sort $OUTDIR/$OUTNAME.bam $OUTDIR/${OUTNAME}.sort
  # mv $OUTDIR/${OUTNAME}.sort.bam $OUTDIR/$OUTNAME.bam
  $samtools index $OUTDIR/${OUTNAME}.sort.bam $OUTDIR/${OUTNAME}.sort.bai
  # rm $OUTDIR/$OUTNAME.sam

  # get idx stats
  echo "$REF_IDX_NAME" >>$OUTDIR/$OUTNAME.info.txt
  $samtools idxstats $OUTDIR/$OUTNAME.sort.bam >>$OUTDIR/$OUTNAME.info.txt
}


# iterate through fasta reference files
# EXCLUDE="human|mouse"
for file in $( ls $REFPATH | grep -Ev "$EXCLUDE"); do
  REFPATHFULL=$REFPATH/$file
  REF_IDX_NAME=$( basename $REFPATHFULL | sed 's/\.f.*$//' )
  # >$OUTDIR/$OUTNAME.info.txt
  map_bwa $REFPATHFULL $REF_IDX_NAME
done

# iterate through EXREF
if [[ -d $EXREFPATH ]]; then
  for file in $( ls $EXREFPATH ); do
    EXREFPATHFULL=$EXREFPATH/$file
    REF_IDX_NAME=$( basename $EXREFPATHFULL | sed 's/\.f.*$//' )
    # >$OUTDIR/$OUTNAME.info.txt
    map_bwa $EXREFPATHFULL $REF_IDX_NAME
  done
fi

# extract contig id/coverage from mapped files and convert to tab seperated table
echo "searchRef_Extract mapping percentage and convert to list"
>$OUTDIR/$OUTNAME.info.list.txt
# post process info file to get tab seperated output
REFNAME_FULL=NA
SUM_GI=NA
SUM_LEN=NA
SUM_READSM=NA
SUM_READSU=NA
SUM_CHROM=-1
cat  $OUTDIR/$OUTNAME.info.txt | while read -r line; do
  LL=$( echo -e "$line" | awk '{ print NF }' )
  if [ "$LL" -eq 1 ]; then 
    REFNAME_FULL=$line
  fi
  if [ "$LL" -gt 1 ]; then
    # set ids
    GI=$( echo -e "$line" | cut -f1 | cut -d"|" -f4 )
    if [ "$SUM_GI" == "NA" ]; then
      SUM_GI=$GI
    else
      SUM_GI=$( echo "$SUM_GI,$GI" )
    fi
    # set length
    LEN=$( echo -e "$line" | cut -f2 )
    if [ "$SUM_LEN" == "NA" ]; then
     SUM_LEN=$LEN
    else
      SUM_LEN=$(( $SUM_LEN + $LEN ))
    fi
    # set mapped
    MAPP=$( echo -e "$line" | cut -f3 )
    if [ "$SUM_READSM" == "NA" ]; then
     SUM_READSM=$MAPP
    else
      SUM_READSM=$(( $SUM_READSM + $MAPP ))
    fi
    # set unmapped
    MAPU=$( echo -e "$line" | cut -f4 )
    if [ "$SUM_READSU" == "NA" ]; then
     SUM_READSU=$MAPU
    else
      SUM_READSU=$(( $SUM_READSU + $MAPU ))
    fi
    # count chromosomes
    SUM_CHROM=$(( $SUM_CHROM+1 ))

  fi
  if [[ "$line" == \** ]]; then
    SUM_GI=$( echo $SUM_GI | sed 's/,.$//' )
    PCT=$( echo -e "$SUM_READSM\t$SUM_READSU" | awk '{ pct=$1*100/($1 + $2); print pct }' )
    echo -e "$REFNAME_FULL\t$SUM_GI\t$SUM_CHROM\t$SUM_LEN\t$SUM_READSM\t$SUM_READSU\t$PCT" >>$OUTDIR/$OUTNAME.info.list.txt
    # reset all
    REFNAME_FULL=NA
    SUM_GI=NA
    SUM_LEN=NA
    SUM_READSM=NA
    SUM_READSU=NA
    SUM_CHROM=-1
  fi
done

# sort according to most mapped
sort -k7,7 -nr $OUTDIR/$OUTNAME.info.list.txt >$OUTDIR/$OUTNAME.info.list.sort.txt

# filter out unwanted refereces for target reference search if BLACKLIST not empty, else just rename file
if [[ ! -z $BLACKLIST ]]; then
  cat $OUTDIR/$OUTNAME.info.list.sort.txt | grep -Ev "$( echo "$BLACKLIST" | tr "\n" "|" | sed 's/|$//' )" >$OUTDIR/$OUTNAME.info.list.sort.filter.txt
else
  cp $OUTDIR/$OUTNAME.info.list.sort.txt $OUTDIR/$OUTNAME.info.list.sort.filter.txt
fi

# write best scoring reference to variant.config
REF_BEST=$( head -n1 $OUTDIR/$OUTNAME.info.list.sort.filter.txt | awk '{ print $1".fasta" }' )
REF_BEST_FULL=${REFPATH}/${REF_BEST}
sed -i 's#^REF_FA_ORIGINAL=.*#REF_FA_ORIGINAL='"${REF_BEST_FULL}"'#' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
sed -i 's#^REF_MAPPING=.*#REF_MAPPING='"${OUTDIR}/exref"'#' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config

# summary
echo "Top_five_best_matching_organisms:" >$OUTDIR/$OUTNAME.info.summary.txt
cat $OUTDIR/$OUTNAME.info.list.sort.txt | head -n5 | cut -f1,3- | awk '{ print "sumref\t"$0 }' >>$OUTDIR/$OUTNAME.info.summary.txt
echo "Best_matching_reference_organism:" >>$OUTDIR/$OUTNAME.info.summary.txt
cat $OUTDIR/$OUTNAME.info.list.sort.filter.txt | head -n1 | cut -f1,3- | awk '{ print "topref\t"$0 }' >>$OUTDIR/$OUTNAME.info.summary.txt
echo "Top_matching_contamination_organism( >=5%):" >>$OUTDIR/$OUTNAME.info.summary.txt
cat $OUTDIR/$OUTNAME.info.list.sort.txt | grep -E $(echo "$BLACKLIST" | tr "\n" "|" | sed 's/|$//') | awk '{if ($7 >= 5) print $0}' | cut -f1,3- | awk '{ print "extref\t"$0 }' >>$OUTDIR/$OUTNAME.info.summary.txt

# write path to alternate/contamination references to variant config and create path with symbolic links
BLACKREF=$( cat "$OUTDIR/$OUTNAME.info.summary.txt" | grep -e 'extref' | cut -f2 )
mkdir -p $OUTDIR/exref
echo "$BLACKREF" | while read -r line; do
  echo $line
  REXNAME=$( ls "$EXREFPATH" | grep -e "$line" )
  REXFULL="$EXREFPATH/$REXNAME"
  ln -s $REXFULL $OUTDIR/exref/$REXNAME
done

# remove temp files
rm $OUTDIR/$OUTNAME*sam
rm $OUTDIR/$OUTNAME.info.list.txt
rm $OUTDIR/$OUTNAME*fastq
rm $OUTDIR/$OUTNAME*bam $OUTDIR/$OUTNAME*bai
