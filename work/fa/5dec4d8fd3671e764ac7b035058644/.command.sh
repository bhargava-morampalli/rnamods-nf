#!/bin/bash -ue
samtools view -S -b -h ivt_3_23s.sam | samtools sort -o ivt_3_23s_sorted.bam
