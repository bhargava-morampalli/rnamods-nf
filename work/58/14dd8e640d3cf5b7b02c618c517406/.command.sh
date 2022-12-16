#!/bin/bash -ue
seqkit fq2fa ivt_1_23s.fastq -o ivt_1_23s.fasta
seqkit seq ivt_1_23s.fasta -n -i > ivt_1_23s.txt
