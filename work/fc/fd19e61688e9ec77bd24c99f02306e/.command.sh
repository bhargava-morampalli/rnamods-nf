#!/bin/bash -ue
samtools view -S -b -h ivt_3_16s.sam | samtools sort -o ivt_3_16s_sorted.bam
