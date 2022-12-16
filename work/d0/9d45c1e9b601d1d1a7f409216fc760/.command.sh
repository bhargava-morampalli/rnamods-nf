#!/bin/bash -ue
minimap2 -ax splice -uf -k14 --secondary=no k12_16s_88_extended.fa ivt_3_16s.fastq | samtools view -S -b -h | samtools sort -o ivt_3_16s_sorted.bam
