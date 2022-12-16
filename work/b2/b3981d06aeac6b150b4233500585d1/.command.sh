#!/bin/bash -ue
samtools view -S -b -h ivt_1_23s.sam | samtools sort -o ivt_1_23s_sorted.bam
