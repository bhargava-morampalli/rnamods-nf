#!/bin/bash -ue
samtools view -S -b -h native_2_16s.sam | samtools sort -o native_2_16s_sorted.bam
