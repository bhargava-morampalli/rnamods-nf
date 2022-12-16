#!/bin/bash -ue
/home/bhargavam/f5c-v0.7/f5c eventalign --scale-events --signal-index --print-read-names --rna -r ivt_1_16s.fastq -b ivt_1_16s_sorted.bam -g k12_16s_88_extended.fa --summary ivt_1_16s_summary.txt --threads 40 > ivt_1_16s_eventalign.txt
