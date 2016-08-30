#!/bin/bash

REGEX_FILTER=$1
OUTNAME=collected_ccf.ccf

# 1. first concatenate all ccf files
>$OUTNAME
for file in $( ls -d */ | grep -Ev 'batch|reference|bwa_index' | grep -E "${REGEX_FILTER}" ); do
  echo "$file" 
  awk '{print}' < $file/vcfs/*.ccf >>$OUTNAME

done

# replace NC_005861_PAC_03_d38_d63_d20_i46_s with NC_005861.PAC03 as reference
sed -i 's|NC_005861_PAC_03_d38_d63_d20_i46_s|NC_005861\.PAC03|g' $OUTNAME

# filter out diverse calls (phiX and PH_S)
>filter_$OUTNAME
cat $OUTNAME | grep -Ev 'phiX_NODE|PH_S' >>filter_$OUTNAME


# 2. use phpmyadmin to import ccf file
# using Format: CSV using DATA LOAD
# Columns separated with: \t
# Columns enclosed with:
# Columns escaped with: \
# Lines terminated with: auto
# Column names: project_id,experiment_id,sequencing_id,experiment,replicate,selection,selection_value,timepoint,filename,iteration,chrom,position,ident,ref,alt,qual,filter,rep_pos,svtype,length,caller,hp_len,hp_seq,cov_total,cov_ref,cov_alt,percentage_alt,snpeff,snpeff_gene

# Delete rows with the following experiment ids
# DELETE FROM `vcf2ccf5` WHERE `experiment_id` = 36931 
