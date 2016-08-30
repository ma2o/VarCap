#!/bin/bash
#$-q all.q@cube[ab]*

# removes the illumina multiplex adapter from the reads, per default the multiplex identifier is set to NNNNNN
# however you can provide the exact sequence within the command line

# input variables
FASTQ=$1
PROJ_NAME=$2
OUT_NAME=$3
MPLEX_UNIQUE_ID=$4

# fixed variables (for multiplexed)
ERROR_RATE=0.02
# illumina trueseq indexed adapter sequence for multiplexing
MPLEX_ADAPTER_SEQ1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
MPLEX_ADAPTER_UNIQUE=NNNNNN
if [[ -n $MPLEX_UNIQUE_ID ]];
then
  MPLEX_ADAPTER_UNIQUE=$MPLEX_UNIQUE_ID
fi
MPLEX_ADAPTER_SEQ2=ATCTCGTATGCCGTCTTCTGCTTG
INDEX_ADAPTER_SEQ=${MPLEX_ADAPTER_SEQ1}${MPLEX_ADAPTER_UNIQUE}${MPLEX_ADAPTER_SEQ2}
echo "Used illumina multiplex indexed adapter:"
echo "$INDEX_ADAPTER_SEQ"
# illumina truseq universal adapter sequence
UNIVERSAL_ADAPTER=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
UNIVERSAL_ADAPTER_REV_COMPL=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
echo "Universal adaptor reverse complement:"
echo "$UNIVERSAL_ADAPTER_REV_COMPL"

# set minimum length
MIN_LENGTH=40

cd ${PROJ_NAME}/filter

### use cutadapt to remove illumina adapers

cutadapt -e $ERROR_RATE -m $MIN_LENGTH -b $INDEX_ADAPTER_SEQ -b $UNIVERSAL_ADAPTER_REV_COMPL $FASTQ >${OUT_NAME}.fastq

