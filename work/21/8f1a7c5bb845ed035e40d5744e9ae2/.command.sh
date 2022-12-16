#!/bin/bash -ue
seqkit fq2fa native_2_23s.fastq -o native_2_23s.fasta
seqkit seq native_2_23s.fasta -n -i > native_2_23s.txt
