#!/bin/bash

IN1=$1
IN2=$2
INLIST=$3
OUT1=exreads_1.fastq.gz
OUT2=exreads_2.fastq.gz

# use filterbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string> include=<t/f>
# from the bbmap package

filterbyname.sh in=$IN1 in2=$IN2 out=$OUT1 out2=$OUT2 names=$INLIST include=t ow=t substring=name truncate=f
