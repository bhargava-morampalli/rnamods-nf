#!/bin/bash -ue
minimap2 -ax splice -uf -k14 --secondary=no k12_23s_78_extended.fa native_3_23s.fastq | samtools view -S -b -h | samtools sort -o native_3_23s_sorted.bam
