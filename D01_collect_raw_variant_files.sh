#!/bin/bash
#$-q all.q@cube[ab]*

CURRENT_DIR=$(pwd)
# open project specific config file
CONFIG_FILE=$CURRENT_DIR/variant.config
. $CONFIG_FILE
# open system/program settings config file in varcap main dir
SYSTEM_CONFIG=${PATH_VARCAP}/program.config
. $SYSTEM_CONFIG

# copy raw calling files to raw_variants folder
RAW_VARIANTS_FOLDER=$PATH_PROJECTS_DATA/$PROJ_NAME/vcfs_raw/raw_variant_files/$BAM_NAME_BASE
mkdir -p $RAW_VARIANTS_FOLDER
# copy gatk files
for file in $PATH_GATK_DATA/$PATH_DATA/$BAM_NAME_BASE*filtered.vcf; 
do 
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/gatk_v${ITER}_$(basename $file) &> /dev/null || echo "gatk variant files not found.";
done
# copy samtools/bcftools files
for file in $PATH_SAMTOOLS_DATA/$PATH_DATA/$BAM_NAME_BASE*.vcf;
do 
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/samtools_v${ITER}_$(basename $file) &> /dev/null || echo "samtools variant files not found.";
done
# copy varscan
for file in $PATH_VARSCAN_DATA/$PATH_DATA/$BAM_NAME_BASE*filter.{snp,indel};
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/varscan_v${ITER}_$(basename $file);
done
# copy lofreq
for file in $PATH_LOFREQ_DATA/$PATH_DATA/$BAM_NAME_BASE*_lofreq_filter;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/lofreq_v${ITER}_$(basename $file);
done
# copy lofreq2
for file in $PATH_LOFREQ2_DATA/$PATH_DATA/$BAM_NAME_BASE*_filter.vcf;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/lofreq2_v${ITER}_$(basename $file);
done
# copy lofreq21
for file in $PATH_LOFREQ21_DATA/$PATH_DATA/$BAM_NAME_BASE*.vcf;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/lofreq21_v${ITER}_$(basename $file);
done
# copy freebayes
for file in $PATH_FREEBAYES_DATA/$PATH_DATA/$BAM_NAME_BASE*_pool.vcf;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/freebayes_v${ITER}_$(basename $file);
done
# copy breakdancer
for file in $PATH_BREAKD_DATA/$PATH_DATA/$BAM_NAME_BASE*.ctx;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/breakd_v${ITER}_$(basename $file);
done
# copy pindel
for file in $PATH_PINDEL_DATA/$PATH_DATA/$BAM_NAME_BASE*_{BP,D,D.vcf,INV,INV.vcf,LI,SI,SI.vcf,TD,TD.vcf};
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/pindel_v${ITER}_$(basename $file);
done
# copy pindel_025
for file in $PATH_PINDEL_DATA_025/$PATH_DATA/$BAM_NAME_BASE*ALL.vcf;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/pindel025_v${ITER}_$(basename $file);
done
# copy delly 0.0.9
for file in $PATH_DELLY_DATA/$PATH_DATA/$BAM_NAME_BASE*.txt;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/delly_v${ITER}_$(basename $file);
done
# copy delly 0.7.2
for file in $PATH_DELLY_072_DATA/$PATH_DATA/$BAM_NAME_BASE*.vcf;
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/delly_072_v${ITER}_$(basename $file);
done
# copy cortex
for file in $PATH_CORTEX_DATA/$PATH_DATA/$BAM_NAME_BASE*/vcfs/$BAM_NAME_BASE*{all_k.raw.vcf,all_k.decomp.vcf};
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/cortex_v${ITER}_$(basename $file) &> /dev/null || echo "cortex variant files not found.";
done
# copy cortex 10521
for file in $( ls $PATH_CORTEX_DATA_10521/$PATH_DATA/$BAM_NAME_BASE*/vcfs/$BAM_NAME_BASE*{all_k.raw.vcf,all_k.decomp.vcf} );
do
  ITER=$( basename $file | grep -o "bwa_[0-9]*[\._]" | grep -o '[0-9]*' )
  cp $file $RAW_VARIANTS_FOLDER/cortex10521_v${ITER}_$(basename $file) &> /dev/null || echo "cortex10521 variant files not found.";
done

echo "Copy of raw variant files finished."

