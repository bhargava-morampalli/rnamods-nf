#!/bin/bash -ue
/home/bhargavam/f5c-v0.7/f5c eventalign --scale-events --signal-index --print-read-names --rna -r ivt_2_23s.fastq -b ivt_2_23s_sorted.bam -g k12_23s_78_extended.fa --summary ivt_2_23s_summary.txt --threads 40 > ivt_2_23s_eventalign.txt
