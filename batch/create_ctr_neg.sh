#!/bin/bash

HUMAN="/mirror/ncbi/current/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.33_GRCh38.p7/GCF_000001405.33_GRCh38.p7_genomic.fna.gz"
zcat "$HUMAN" >/scratch/zojer/projects/genomes/ch_tra_neg/human.fasta
MOUSE="/mirror/ncbi/current/genomes/refseq/vertebrate_mammalian/Mus_musculus/reference/GCF_000001635.24_GRCm38.p4/GCF_000001635.24_GRCm38.p4_genomic.fna.gz"
zcat "$MOUSE" >/scratch/zojer/projects/genomes/ch_tra_neg/mouse.fasta
