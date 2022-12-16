#!/bin/bash -ue
samtools view -S -b -h native_2_23s.sam | samtools sort -o native_2_23s_sorted.bam
