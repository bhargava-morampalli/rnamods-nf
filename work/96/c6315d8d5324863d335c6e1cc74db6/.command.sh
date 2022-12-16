#!/bin/bash -ue
minimap2 -ax splice -uf -k14 k12_16s_88_extended.fa ivt_3.fastq > ivt_3_16s.sam --secondary=no
