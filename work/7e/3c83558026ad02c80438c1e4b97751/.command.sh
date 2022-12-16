#!/bin/bash -ue
samtools view -S -b -h native_3_16s.sam | samtools sort -o native_3_16s_sorted.bam
