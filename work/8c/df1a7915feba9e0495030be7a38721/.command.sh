#!/bin/bash -ue
samtools view -S -b -h native_1_16s.sam | samtools sort -o native_1_16s_sorted.bam
