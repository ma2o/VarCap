#!/usr/bin/python

from Bio import SeqIO
import sys

#queryFile = sys.argv[1]
#try:
#  suppressQualFile = sys.argv[2]
#except IndexError:
#  suppressQualFile = False

#print >> sys.stderr, "converting %s..." % queryFile

#if queryFile.endswith('.fastq'):
#  fastafile = queryFile[:-len('.fastq')] + ".fasta"
#  qualfile = queryFile[:-len('.fastq')] + ".qual"
#elif queryFile.endswith('.fq'):
#  fastafile = queryFile[:-len('.fq')] + ".fasta"
#  qualfile = queryFile[:-len('.fq')] + ".qual"
#else:
#  fastafile = queryFile + ".fasta"
#  qualfile = queryFile + ".qual"

#SeqIO.convert(queryFile, "fastq", fastafile, "fasta")
SeqIO.convert(sys.stdin, "fastq", sys.stdout, "fasta")

#if not suppressQualFile:
#  SeqIO.convert(queryFile, "fastq", qualfile, "qual")
  #SeqIO.convert(sys.stdin, "fastq", sys.stdout, "qual")

