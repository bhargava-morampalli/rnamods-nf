#!/bin/bash -ue
samtools view -S -b -h ivt_2_16s.sam | samtools sort -o ivt_2_16s_sorted.bam
