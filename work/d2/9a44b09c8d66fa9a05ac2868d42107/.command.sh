#!/bin/bash -ue
seqkit fq2fa ivt_3_23s.fastq -o ivt_3_23s.fasta
seqkit seq ivt_3_23s.fasta -n -i > ivt_3_23s.txt
