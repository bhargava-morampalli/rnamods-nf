#!/bin/bash -ue
seqkit fq2fa native_2_16s.fastq -o native_2_16s.fasta
seqkit seq native_2_16s.fasta -n -i > native_2_16s.txt
