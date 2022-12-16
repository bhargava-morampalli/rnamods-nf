#!/bin/bash -ue
minimap2 -ax splice -uf -k14 k12_16s_88_extended.fa ivt_2.fastq > ivt_2_16s.sam --secondary=no
