#!/bin/bash -ue
samtools view -S -b -h ivt_2_23s.sam | samtools sort -o ivt_2_23s_sorted.bam
