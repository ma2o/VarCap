#!/bin/bash

# SLURM
#SBATCH --job-name=wVCbwamem
#SBATCH --cpus-per-task=2
#SBATCH --mem=10000
#SBATCH --nice=1000
#SBATCH --output=bwamem2bam2-%A_%a.out
#SBATCH --error=bwamem2bam2-%A_%a.err

# . /etc/profile

# run bwa sam pipeline
PROJ_BASE=$1
IT=$2
CONFIG_FILE=$PROJ_BASE/variant.config
. $CONFIG_FILE
SYSTEM_CONFIG=$3
. $SYSTEM_CONFIG

TEMP=$TMPDIR
# TEMP=/scratch/zojer/tmpdir

OUTFILE=${BAM_NAME_BASE}_bwa_${IT}
REF=${REF_FA}
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
# BS=$( cat $PATH_PROJECTS_DATA/${PROJ_NAME}/info.txt | grep -e 'BothS' | cut -d' ' -f2 )
BS=$( zcat $PATH_ALT_READS1 | awk 'END {count=NR/4; print count }' )
if [ "$SUBSAMPLE_SIZE_ALT" -lt 0 ]; then
    sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$BS"'/' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
    echo "VARCAP_WARN: No readcount defined, $BS available, taking all reads available."
    echo "VARCAP_WARN: No readcount defined, $BS available, taking all reads available." >>$PATH_LOGS/log.txt
    SUBSAMPLE_SIZE_ALT=$BS
elif [ "$SUBSAMPLE_SIZE_ALT" -gt "$BS" ]; then
    sed -i 's/SUBSAMPLE_SIZE_ALT=.*/SUBSAMPLE_SIZE_ALT\='"$BS"'/' $PATH_PROJECTS_DATA/${PROJ_NAME}/variant.config
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available."
    echo "VARCAP_WARN: Need $SUBSAMPLE_SIZE_ALT reads, but only $BS available, taking all reads available." >>$PATH_LOGS/log.txt
    SUBSAMPLE_SIZE_ALT=$BS
fi

REFNAME=$(basename $REF_FA | sed 's/\.f.*$//')
# subsample simulated reads, if not all reads are needed (path to reference/variant reads set in subscript)
if [ "$SUBSAMPLE_SIZE_ALT" -lt "$BS" ]; then
  mkdir -p $PATH_SUBSAMPLE_DATA
  cd $PATH_SUBSAMPLE_DATA
  bash ${PATH_SCRIPTS}/22_sample_subsets.sh $SUBSAMPLE_SIZE_ALT $BAM_NAME_BASE $REPEATS $PATH_ALT_READS1 $PATH_ALT_READS2 $PATH_SCRIPTS $IT
else
  # make symbolic links instead of subsample output
  mkdir -p $PATH_SUBSAMPLE_DATA
  rm $PATH_SUBSAMPLE_DATA/*.fq
  ln -s ${PATH_ALT_READS1} ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${IT}_1.fq
  ln -s ${PATH_ALT_READS2} ${PATH_SUBSAMPLE_DATA}/${BAM_NAME_BASE}_r${IT}_2.fq
fi

##### 2. Functions run bwa/indexing on EXREF(contamination references) first and take the unmapped reads for mapping to desired reference

ref_index(){
  local REFNAME=$1
  local REFPATH=$2
  ### 2.1 check if index exists, else build it
  if [ -d "${PATH_BWA_075_INDEX}/bwa_index_${REFNAME}" ]; then
      echo "VARCAP:bwa_index exists, using it."
  else
      echo "VARCAP:constructing bwa index."
      mkdir -p $PATH_BWA_075_INDEX/bwa_index_${REFNAME}
      cd $PATH_BWA_075_INDEX
      bash $PATH_SCRIPTS/A01_bwa_build_index.sh ${REFPATH}/${REFNAME} $PATH_BWA_075 $PATH_BWA_075_INDEX $PATH_LOGS
      # qsub -l vf=5G -pe parallel 2 $PATH_SCRIPTS/A01_bwa_build_index.sh ${REF_FA}_map.fasta $PATH_BWA_075 $PATH_BWA_075_INDEX
  fi
}

map_bwa() {
	local REFNAME=$1
	local REFPATH=$2
	local OUTFILE=$3
	local READ1=$4
	local READ2=$5
	
	# cd $PATH_BWA_DATA
	echo "VARCAP: align with bwa mem"
        echo "VARCAP: align with bwa mem" >>$PATH_LOGS/logs.txt
        	
	# copy index to $TEMP dir
	mkdir -p $TEMP/bwa_index
	cp -r $PATH_BWA_075_INDEX/bwa_index_$REFNAME  $TEMP/bwa_index/
	
	# check if READ1 input is bam, else fastq
	if [[ "$READ1" == *.bam ]]; then
	  # bam
	  $PATH_SAMTOOLS/samtools bam2fq $READ1 | ${PATH_BWA_075}/bwa mem -M -t 2 -p $TEMP/bwa_index/bwa_index_$REFNAME/$REFNAME - >$TEMP/$OUTFILE.sam
	else
	  # second align reads M=mark secondary,p=paired, -T 10=min qual score
          echo "VARCAP: bwa mem properties:" >>$PATH_LOGS/logs.txt
          echo "VARCAP: $PATH_BWA_075/bwa mem -M -t 2 -p $TEMP/bwa_index/bwa_index_$REFNAME/$REFNAME $READ1 $READ2" >>$PATH_LOGS/logs.txt
	  $PATH_BWA_075/bwa mem -M -t 2 $TEMP/bwa_index/bwa_index_$REFNAME/$REFNAME $READ1 $READ2 >$TEMP/$OUTFILE.sam
	fi
	
	# sam to bam conversion and sorting
	$PATH_SAMTOOLS/samtools view -bS $TEMP/$OUTFILE.sam | $PATH_SAMTOOLS/samtools sort - $TEMP/$OUTFILE
	
	# indexing bam file
	$PATH_SAMTOOLS/samtools index $TEMP/$OUTFILE.bam $TEMP/$OUTFILE.bai
	
	# removing .sam and .sai files
	if [ -f $TEMP/$OUTFILE.bai ];
	  then
	  rm $TEMP/$OUTFILE.sam $TEMP/$READNAME1.sai $TEMP/$READNAME2.sai
	  rm $TEMP/run_bwa2bam.sh.po* $TEMP/run_bwa2bam.sh.pe*
	fi
}

##### 3. Main 

### 3.1 Check index and map
# iterate through contamination references
if [ -d "$REF_MAPPING" ]; then
  # echo "Iterating contamination REF_MAPPING: $EXNAME1"
  R1="$READS1"
  R2="$READS2"
  for EXNAME in $( ls "$REF_MAPPING" ); do
    EXNAME1=$( basename "$EXNAME" | sed 's/\.f.*a$//' )
    echo "Iterating contamination REF_MAPPING: $EXNAME1"
    ref_index $EXNAME1 $REF_MAPPING
    map_bwa $EXNAME1 $REF_MAPPING $EXNAME1 $R1 $R2
    # extract unmapped reads
    $PATH_SAMTOOLS/samtools view -b -f 4 ${TEMP}/${EXNAME1}.bam >${TEMP}/${EXNAME1}_unmapped.bam
    mv ${TEMP}/${EXNAME1}_unmapped.bam ${TEMP}/unmapped.bam
    R1="${TEMP}/unmapped.bam"
    R2="NONE"
    rm ${TEMP}/${EXNAME1}.bam
    cp ${TEMP}/unmapped.bam ${PATH_BWA_DATA}
  done
fi

# check/index target reference
if [ -f "$REF" ]; then
  # echo "Processing target REF: $REFNAME"
  RT1="$READ1"
  RT2="$READ2"
  # check if unmapped bam exists
  if [ -f "${PATH_BWA_DATA}/unmapped.bam" ]; then
    RT1="${PATH_BWA_DATA}/unmapped.bam"
    RT2="NONE"  
  fi
  REFNAME=$( basename $REF | sed 's/\.f.*a$//' )
  REFDIR=$( dirname $REF | sed 's/\/$//' )
  echo "Processing target REF: $REFNAME"
  echo "Check index:"
  ref_index $REFNAME $REFDIR
  echo "Map with bwa:"
  map_bwa $REFNAME $REFDIR $OUTFILE $RT1 $RT2
  # $PATH_SAMTOOLS/samtools index $TEMP/$OUTFILE.bam $TEMP/$OUTFILE.bai
  cp ${TEMP}/$OUTFILE.bam ${TEMP}/$OUTFILE.bai ${PATH_BWA_DATA}
fi


### calculate parameters for mapping file (insert size,coverage)
# print insert size
java -jar $PATH_PICARD/picard.jar CollectInsertSizeMetrics I=$PATH_BWA_DATA/$OUTFILE.bam H=$PATH_BWA_DATA/${OUTFILE}.IS.pdf O=$PATH_BWA_DATA/${OUTFILE}.IS.txt
MAP_STATS=$($PATH_SAMTOOLS/samtools depth $PATH_BWA_DATA/$OUTFILE.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "map_average_cov\t",sum/NR; print "map_stdev_cov\t",sqrt(sumsq/NR - (sum/NR)**2)}')
echo "$MAP_STATS" >$PATH_BWA_DATA/${OUTFILE}.COV.txt
echo "$MAP_STATS" >>$PATH_PROJECTS_DATA/${PROJ_NAME}/info.txt

# log end
echo "$JOB_ID finished" >$PATH_LOGS/D02_${SLURM_JOBID}.txt
