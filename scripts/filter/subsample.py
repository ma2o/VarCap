#!/usr/bin/env python
import argparse
import random
import subprocess
import sys

parser = argparse.ArgumentParser(description=\
  'This program outputs to stdout a user-defined number of reads from given reads file[s]')
parser.add_argument('-pe',action="store_true",dest="pairedEnd",required=False,\
  help="specifies if paired-end (default is single-end)")
parser.add_argument('-fasta',action="store_true",dest="isFasta",default=False,required=False,\
  help="specifies if fasta (default is fastq)")
parser.add_argument('-n',action="store",type=int,dest="numSeq",required=True,\
  help="specifies number of reads to output")
parser.add_argument('-in',action="store",dest="inputs",nargs='+',required=True,\
  help="specifies input files (put both input files separated by a space for paired end)")
parser.add_argument('-out',action="store",dest="outputs",nargs='*',required=False,\
  default=['/dev/stdout'],help="specifies output files (default for single end is stdout)")
args = parser.parse_args()
isFasta = args.isFasta
numSeq = args.numSeq
pairedEnd = args.pairedEnd
reads = args.inputs
outs = args.outputs
if not(pairedEnd):
  if len(reads) != 1 or len(outs) > 1:
    print "please specify only one read file and at most one output file (if paired end, specify with -pe)"
    sys.exit()
else:
  if len(reads) != 2 or len(outs) != 2:
    print "please specify both read files and two output files"
    sys.exit()
fileName = reads[0]
if(fileName.find(".gz") != -1):
  subprocess.call("gunzip "+fileName, shell=True)
  fileName = fileName[0: -3]
random.seed()
#name of the input file (fasta or fastq)
#assumes input file is standard fasta/fastq format
increment = 0

#if it's a fasta file
if isFasta:
  increment = 2
#else if it's a fastq file
else:
  increment = 4
ttl = sum(1 for line in open(fileName))/increment
if(ttl < numSeq):
  sys.stdout.write("You only have "+str(ttl)+" reads!\n")
  sys.exit()
totalReads = range(0, ttl)
random.shuffle(totalReads)
totalReads = totalReads[0: numSeq]
totalReads.sort()

FILE = open(fileName, 'r')
listIndex = 0

if not(pairedEnd):
  out = open(outs[0], 'w')
  for i in range(0, ttl):
    curRead = ""
    for j in range(0, increment):
      curRead += FILE.readline()
    if (i == totalReads[listIndex]):
      out.write(curRead)
      listIndex += 1
      if(listIndex == numSeq):
        break
  FILE.close()
  out.close()
else:
  fileName2 = reads[1]
  if(fileName2.find(".gz") != -1):
    subprocess.call("gunzip "+fileName2, shell=True)
    fileName2 = fileName2[0: -3]
  out = open(outs[0], 'w')
  out2 = open(outs[1], 'w')
  FILE2 = open(fileName2, 'r')
  for i in range(0, ttl):
    curRead = ""
    curRead2 = ""
    for j in range(0, increment):
      curRead += FILE.readline()
      curRead2 += FILE2.readline()
    if (i == totalReads[listIndex]):
      out.write(curRead)
      out2.write(curRead2)
      listIndex += 1
      if(listIndex == numSeq):
        break
  FILE.close()
  FILE2.close()
  out.close()
  out2.close()