#!/usr/bin/env python

import sys

namesfilename = sys.argv[1]
is_blacklist  = int(sys.argv[2])

def write_fastq(lines):
  name = lines[0][1:-1].strip()
  if name.endswith("/1") or name.endswith("/2"):
    name = name[:-2]
  if is_blacklist and not names.has_key(name) or not is_blacklist and names.has_key(name):
    for line in lines:
      sys.stdout.write(line)

names={}
infile=open(namesfilename)
for line in infile:
  names[line.strip()]=1
infile.close()

lines=[]
for line in sys.stdin:
  lines.append(line)
  if len(lines) == 4:
    write_fastq(lines)
    lines=[]
for line in lines:
  sys.stderr.write(line)
