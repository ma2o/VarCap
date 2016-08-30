#!/bin/bash
#$-q all.q@cube[ab]*

# INPUT_FILES=$@;
INPUT_FILE_DIR=$1

INPUT_FILES=${INPUT_FILE_DIR}/*.ccf
PASS=$2

DB_HOST=vdbdev
DB_DATABASE=varcap_test
DB_USER=varcap_oo1
DB_TABLE=vcf2ccf3

COUNTER=1
for file in $INPUT_FILES;do
  echo "insert file: "$COUNTER
  mysql --host=$DB_HOST --user=$DB_USER --password=$PASS $DB_DATABASE -e "LOAD DATA LOCAL INFILE '${file}' INTO TABLE ${DB_TABLE} ( project_id,experiment_id,sequencing_id,experiment,replicate,selection,selection_value,timepoint,filename,iteration,chrom,position,ident,ref,alt,qual,filter,rep_pos,svtype,length,caller,hp_len,hp_seq,cov_total,cov_ref,cov_alt,percentage_alt,snpeff,snpeff_gene);
  SELECT COUNT(*) FROM ${DB_TABLE};"
  let COUNTER=COUNTER+1
done
PASS=none



