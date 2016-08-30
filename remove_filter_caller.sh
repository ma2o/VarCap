#!bin/bash
#$-q all.q@cube[ab]*

FOLDER_NAME=$1
BATCH_COMMAND=$2
COMMAND_PROP=$3

CURR_DIR=$( pwd )

if [ "$#" -lt 1 ];then
  echo "Usage: bash patch_*.sh <folder_name_substring> <command>"
  exit
fi

# for folder in $( ls | grep -e 'evochlamy_01' ); do rm $folder/vcfs_raw/*.vcf; rm $folder/vcfs_raw/vcfs_temp/*.vcf; done

for folder in $( ls | grep -e "$FOLDER_NAME" );
do
  # echo $folder
  rm $folder/filter/*.fastq
  rm $folder/mapper/subsample/*
  rm $folder/mapper/subsample/subsets/*
  rm $folder/caller/breakdancer/data/*
  rm -r $folder/caller/cortex/data/*
  rm $folder/caller/delly/data/*
  rm $folder/caller/gatk/data/*
  rm $folder/caller/lofreq/data/*
  rm $folder/caller/lofreq2/data/*
  rm $folder/caller/pindel/data/*
  rm $folder/caller/pindel_025/data/*
  rm $folder/caller/samtools/data/*
  rm $folder/caller/varscan/data/*
  cd $CURR_DIR
done

