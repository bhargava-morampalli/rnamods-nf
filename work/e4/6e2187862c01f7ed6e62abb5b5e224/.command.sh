#!/bin/bash -ue
minimap2 -ax splice -uf -k14 --secondary=no k12_23s_78_extended.fa ivt_2_23s.fastq | samtools view -S -b -h | samtools sort -o ivt_2_23s_sorted.bam
