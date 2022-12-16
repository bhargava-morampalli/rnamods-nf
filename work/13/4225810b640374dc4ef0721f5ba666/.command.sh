#!/bin/bash -ue
seqkit fq2fa native_1_16s.fastq -o native_1_16s.fasta
seqkit seq native_1_16s.fasta -n -i > native_1_16s.txt
