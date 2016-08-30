#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCvs2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=varscan-%A_%a.out
#SBATCH --error=varscan-%A_%a.err

source /etc/profile

REF=$1
BAM=$2
OUTDIR=$3
BAM_BASE=$(basename $BAM | sed 's/\..am//')

PATH_VARSCAN=$4
PATH_SAMTOOLS=$5
#PATH_VARSCAN_DATA=/scratch/zojer/projects/test_pipeline/mapping/snp_calling/varscan

cd $OUTDIR
#sort bam file
$PATH_SAMTOOLS/samtools sort $BAM $BAM_BASE.sorted 

#create pileup file
$PATH_SAMTOOLS/samtools mpileup -C50 -d1000 -Ef $REF $BAM_BASE.sorted.bam >$BAM_BASE.mpileup

#run pileup2snp
java -Xmx8g -jar $PATH_VARSCAN pileup2snp $BAM_BASE.mpileup -min-var-freq 0.001 >$BAM_BASE.snp

#run pileup2indel
java -Xmx8g -jar $PATH_VARSCAN pileup2indel $BAM_BASE.mpileup -min-var-freq 0.001 >$BAM_BASE.indel

#run somatic
#java -Xmx8G -jar $PATH_VARSCAN somatic $BAM_BASE.mpileup $BAM_BASE.somatic --mpileup 1

#filter --min-strands2 1
java -jar $PATH_VARSCAN filter $BAM_BASE.indel --min-avg-qual 25 --min-var-freq 0.001  --output-file $BAM_BASE.filter.indel
java -jar $PATH_VARSCAN filter $BAM_BASE.snp --min-avg-qual 25 --min-var-freq 0.001 --min-strands2 2 --indel-file $BAM_BASE.filter.indel --output-file $BAM_BASE.filter.snp

#delete pileup file
rm $BAM_BASE.mpileup $BAM_BASE.sorted.bam
