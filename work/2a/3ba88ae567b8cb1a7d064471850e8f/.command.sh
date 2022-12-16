#!/bin/bash -ue
samtools view -S -b -h native_1_23s.sam | samtools sort -o native_1_23s_sorted.bam
