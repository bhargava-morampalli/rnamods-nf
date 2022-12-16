#!/bin/bash -ue
seqkit fq2fa ivt_1_16s.fastq -o ivt_1_16s.fasta
seqkit seq ivt_1_16s.fasta -n -i > ivt_1_16s.txt
