#!/bin/bash
#$-q all.q@cube[ab]*

# run bwa sam pipeline
PROJ_BASE=$1
IT=$2
CONFIG_FILE=$PROJ_BASE/variant.config
. $CONFIG_FILE
SYSTEM_CONFIG=$3
. $SYSTEM_CONFIG

OUTFILE=${BAM_NAME_BASE}_bwa_${IT}
REF=${REF_FA}_map.fasta
REFNAME=$( basename $REF | sed 's/\.f.*$//')
# output of subsampling/input for bwa
READ1=${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${IT}_1.fq
READ2=${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${IT}_2.fq
READNAME1=$(basename $READ1 | sed 's/\.f.*$//')
READNAME2=$(basename $READ2 | sed 's/\.f.*$//')


# log start
echo "$JOB_ID started" >$PATH_LOGS/D02_${IT}_${JOB_ID}.txt

### 1. subsampling of readcount

# check if we have enough reads for subsampling, else take everything we have
BS=$( cat $PATH_PROJECTS_DATA/${PROJ_NAME}/info.txt | grep -e 'BothS' | cut -d' ' -f2 )
if [ "$SUBSAMPLE_SIZE_ALT" -gt "$BS" ]; then
    sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$BS"'/' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available."
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available." >>$PATH_LOGS/log.txt
    SUBSAMPLE_SIZE_ALT=$BS
fi

REFNAME=$(basename $REF_FA | sed 's/\.f.*$//')
# subsample simulated reads (path to reference/variant reads set in subscript)
mkdir -p $PATH_SUBSAMPLE_DATA
cd $PATH_SUBSAMPLE_DATA
bash ${PATH_SCRIPTS}/22_sample_subsets.sh $SUBSAMPLE_SIZE_ALT $BAM_NAME_BASE $REPEATS $PATH_ALT_READS1 $PATH_ALT_READS2 $PATH_SCRIPTS $IT


### 2. Check if index exists

# check if index exists, else build it
if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${REFNAME}" ]; then
  echo "VARCAP:bwa_index exists, using it."
else
  echo "VARCAP:constructing bwa index."
  mkdir -p $PATH_BWA_075_INDEX/bwa_index_${REFNAME}
  cd $PATH_BWA_075_INDEX
  bash $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA}_map.fasta $PATH_BWA_075 $PATH_BWA_075_INDEX
  # qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA}_map.fasta $PATH_BWA_075 $PATH_BWA_075_INDEX
fi


### 3. run bwa

cd $PATH_BWA_DATA
echo "VARCAP: align with bwa mem"
# second align reads M=mark secondary,p=paired, -T 10=min qual score
$PATH_BWA_075/bwa mem -M -t 2 $PATH_BWA_075_INDEX/bwa_index_$REFNAME/$REFNAME $READ1 $READ2 >$OUTFILE.sam

echo "build samfile"
# third build sam file
# ${PATH_BWA}/bwa sampe -f $PATH_BWA_DATA/$OUTFILE.sam $PATH_BWA_INDEX/bwa_index_$REFNAME/$REFNAME $READNAME1.sai $READNAME2.sai $READ1 $READ2

# sam to bam conversion and sorting
$PATH_SAMTOOLS/samtools view -bS $PATH_BWA_DATA/$OUTFILE.sam | $PATH_SAMTOOLS/samtools sort - $PATH_BWA_DATA/$OUTFILE

# indexing bam file
$PATH_SAMTOOLS/samtools index $PATH_BWA_DATA/$OUTFILE.bam $PATH_BWA_DATA/$OUTFILE.bai

# removing .sam and .sai files
if [ -f $PATH_BWA_DATA/$OUTFILE.bai ];
  then
  rm $OUTFILE.sam $READNAME1.sai $READNAME2.sai
  rm run_bwa2bam.sh.po* run_bwa2bam.sh.pe*
fi

# split bam files according to references ( create bam for ref, but remove those just for mapping)
echo "VARCAP: split bam to match"
cat $REF | grep -e '^>' | sed 's/^>//' >true_contigs.txt
mv ${OUTFILE}.bam ${OUTFILE}.map.bam
mv ${OUTFILE}.bai ${OUTFILE}.map.bai
# touch ${OUTFILE}.bam
while read line; do
  $PATH_SAMTOOLS/samtools view -b ${OUTFILE}.map.bam $line >${OUTFILE}.temp.bam
  
  if [[ -e ${OUTFILE}.bam ]]; then
    $PATH_SAMTOOLS/samtools cat -o ${OUTFILE}.concat.bam ${OUTFILE}.bam ${OUTFILE}.temp.bam
    mv ${OUTFILE}.concat.bam ${OUTFILE}.bam
    rm ${OUTFILE}.temp.bam
  else
    mv ${OUTFILE}.temp.bam ${OUTFILE}.bam
  fi
  # echo "$PATH_SAMTOOLS/samtools view -b ${OUTFILE}.map.bam $line >>${OUTFILE}.bam"
done <true_contigs.txt


# create unmapped bam
# $PATH_SAMTOOLS/samtools view -b -F 4 ${OUTFILE}.map.bam >${OUTFILE}.F4mapped.bam
$PATH_SAMTOOLS/samtools view -b -f 4 ${OUTFILE}.map.bam >${OUTFILE}.unmapped.bam

# index bam file
$PATH_SAMTOOLS/samtools index $PATH_BWA_DATA/$OUTFILE.bam $PATH_BWA_DATA/$OUTFILE.bai


### calculate parameters for mapping file (insert size,coverage)
# print insert size
java -jar $PATH_PICARD/CollectInsertSizeMetrics.jar I=$PATH_BWA_DATA/$OUTFILE.bam H=$PATH_BWA_DATA/${OUTFILE}.IS.pdf O=$PATH_BWA_DATA/${OUTFILE}.IS.txt
MAP_STATS=$($PATH_SAMTOOLS/samtools depth $PATH_BWA_DATA/$OUTFILE.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "map_average_cov\t",sum/NR; print "map_stdev_cov\t",sqrt(sumsq/NR - (sum/NR)**2)}')
echo "$MAP_STATS" >$PATH_BWA_DATA/${OUTFILE}.COV.txt
echo "$MAP_STATS" >>../../info.txt

# log end
echo "$JOB_ID finished" >$PATH_LOGS/D02_${JOB_ID}.txt
