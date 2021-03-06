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
  cd $folder
  if [ -n $COMMAND_PROP ];then
    echo "executing $BATCH_COMMAND $COMMAND_PROP in $folder"
    bash $BATCH_COMMAND $COMMAND_PROP
  else
    echo "executing $BATCH_COMMAND in $folder"
    bash $BATCH_COMMAND
  fi
  cd $CURR_DIR
done

