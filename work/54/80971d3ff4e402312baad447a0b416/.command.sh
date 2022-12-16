#!/bin/bash -ue
seqkit fq2fa native_3_16s.fastq -o native_3_16s.fasta
seqkit seq native_3_16s.fasta -n -i > native_3_16s.txt
