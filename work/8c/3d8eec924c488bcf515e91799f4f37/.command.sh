#!/bin/bash -ue
/home/bhargavam/f5c-v0.7/f5c eventalign --scale-events --signal-index --print-read-names --rna -r native_2_16s.fastq -b native_2_16s_sorted.bam -g k12_16s_88_extended.fa --summary native_2_16s_summary.txt --threads 40 > native_2_16s_eventalign.txt
