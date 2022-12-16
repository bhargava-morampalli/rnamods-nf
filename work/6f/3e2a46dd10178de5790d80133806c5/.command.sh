#!/bin/bash -ue
samtools view -S -b -h ivt_1_16s.sam | samtools sort -o ivt_1_16s_sorted.bam
