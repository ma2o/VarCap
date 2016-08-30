#!/bin/bash

# SLURM
#SBATCH --job-name=wVCCsam
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=samtools-%A_%a.out
#SBATCH --error=samtools-%A_%a.err

. /etc/profile

REF=$1
SAM=$2
SAM_NAME=$(basename $SAM | sed 's/\..am$//')
SAM_END=$(basename $SAM | sed 's/.*\.\(.am\)$/\1/')
PATH_SAMTOOLS=$3


#if input is sam then convert to bam
if [ $SAM_END == "sam" ]
  then
  #convert bam to sam
  #samtools view -h $SAM > $NAME.bam
  #sam to bam
  ${PATH_SAMTOOLS}/samtools view -bS -o $SAM_NAME.bam $SAM
  SAM=$SAM_NAME.bam
  echo $SAM
fi
#run samtools/bcftools pipeline
echo "sort"
${PATH_SAMTOOLS}/samtools sort $SAM $SAM_NAME.sorted
echo "mpileup"
${PATH_SAMTOOLS}/samtools mpileup -C 50 -d 2000 -m 3 -F 0.0002 -Euf $REF $SAM_NAME.sorted.bam | ${PATH_SAMTOOLS}/bcftools/bcftools view -bvcg - >$SAM_NAME.raw.bcf  
${PATH_SAMTOOLS}/bcftools/bcftools view $SAM_NAME.raw.bcf | ${PATH_SAMTOOLS}/bcftools/vcfutils.pl varFilter -d 5 -D 2000 >$SAM_NAME.snp.flt.vcf

# possible tryouts
# ${PATH_SAMTOOLS}/samtools mpileup -C50 -d1000 -m3 -F0.0002 -Euf $REF $SAM_NAME.sorted.bam | ${PATH_SAMTOOLS}/bcftools/bcftools view -bvcg - > $SAM_NAME.raw.bcf
# ${PATH_SAMTOOLS}/bcftools/bcftools view $SAM_NAME.raw.bcf | ${PATH_SAMTOOLS}/bcftools/vcfutils.pl varFilter -D300 > $SAM_NAME.snp.flt.vcf

