#!/bin/bash -ue
seqkit fq2fa ivt_2_16s.fastq -o ivt_2_16s.fasta
seqkit seq ivt_2_16s.fasta -n -i > ivt_2_16s.txt
