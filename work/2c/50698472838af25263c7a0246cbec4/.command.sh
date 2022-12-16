#!/bin/bash -ue
minimap2 -ax splice -uf -k14 k12_23s_78_extended.fa native_2.fastq > native_2_23s.sam --secondary=no
