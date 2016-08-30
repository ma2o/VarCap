/bin/bash

# script for usage of vmatch and mkvtree
# first create an index using mkvtree
# then run vmatch with previously created index

# mkvtree
mkvtree -db ../NC_005861.fasta -v -pl -sti1 -bwt -dna -bck -suf -lcp -tis -ois -skp

# vmatch
vmatch -v -l 250 ../mkvtree/NC_005861.fasta
vmatch -v -l 250 -sort ia -best 40 -absolute ../mkvtree/NC_005861.fasta
# hamming distance: h (differences: 2 substitutions)
vmatch -v -l 250 -h 2 -sort ia -best 40 -identity 95 -absolute ../mkvtree/NC_005861.fasta
# edit distance: e (differences: substitutions, insertions, deletions)
vmatch -v -l 250 -e 2 -sort ia -best 40 -identity 95 -absolute ../mkvtree/NC_005861.fasta

# test case
vmatch -v -l 250 -h 2 -absolute ../mkvtree/NC_005861.fasta

