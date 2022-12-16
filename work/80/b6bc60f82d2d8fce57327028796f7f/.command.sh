#!/bin/bash -ue
seqkit fq2fa native_3_23s.fastq -o native_3_23s.fasta
seqkit seq native_3_23s.fasta -n -i > native_3_23s.txt
