#!/bin/bash
FILE1=$1
FILE2=$2
SCRIPTDIR=$(dirname $0)
echo $SCRIPTDIR

OUTPREFIX1=$(echo $FILE1 | sed "s/.fastq$//")
OUTPREFIX2=$(echo $FILE2 | sed "s/.fastq$//")

#now get all reads from file one
grep '^@.\+/[12]$' $FILE1 | sed "s/^@//" | sed "s/\/[12]$//"  > $OUTPREFIX1\_temp1.list
#grep '^@.\+\s[12]:' $FILE1 | sed "s/^@//" | sed "s/\s[12]//"  > $OUTPREFIX1\_temp1.list

#reduce file2 to temp1.list
cat $FILE2  | python $SCRIPTDIR/filter_fastqfile.py $OUTPREFIX1\_temp1.list 0 > $OUTPREFIX2\_rp.fastq

#now get all read from FILE2_rp
grep '^@.\+/[12]$' $OUTPREFIX2\_rp.fastq | sed "s/^@//" | sed "s/\/[12]$//" > $OUTPREFIX2\_temp2.list
#grep '^@.\+\s[12]:' $OUTPREFIX2\_rp.fastq | sed "s/^@//" | sed "s/\s[12]//" > $OUTPREFIX2\_temp2.list

#reduce file1 to temp2.list
cat $FILE1 | python $SCRIPTDIR/filter_fastqfile.py $OUTPREFIX2\_temp2.list 0 > $OUTPREFIX1\_rp.fastq

#delete temp files
rm $OUTPREFIX1\_temp1.list
rm $OUTPREFIX2\_temp2.list
