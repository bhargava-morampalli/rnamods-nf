#!/bin/bash -ue
minimap2 -ax splice -uf -k14 k12_16s_88_extended.fa native_1.fastq > native_1_16s.sam --secondary=no
