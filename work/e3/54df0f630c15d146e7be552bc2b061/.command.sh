#!/bin/bash -ue
samtools view -S -b -h native_3_23s.sam | samtools sort -o native_3_23s_sorted.bam
