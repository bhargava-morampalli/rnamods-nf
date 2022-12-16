#!/bin/bash -ue
seqkit fq2fa ivt_2_23s.fastq -o ivt_2_23s.fasta
seqkit seq ivt_2_23s.fasta -n -i > ivt_2_23s.txt
