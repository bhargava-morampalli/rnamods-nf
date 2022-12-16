#!/bin/bash -ue
seqkit fq2fa native_1_23s.fastq -o native_1_23s.fasta
seqkit seq native_1_23s.fasta -n -i > native_1_23s.txt
