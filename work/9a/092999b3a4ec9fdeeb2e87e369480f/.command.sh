#!/bin/bash -ue
seqkit fq2fa ivt_3_16s.fastq -o ivt_3_16s.fasta
seqkit seq ivt_3_16s.fasta -n -i > ivt_3_16s.txt
