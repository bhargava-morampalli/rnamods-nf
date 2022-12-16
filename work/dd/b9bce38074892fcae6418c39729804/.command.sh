#!/bin/bash -ue
tombo detect_modifications level_sample_compare     --fast5-basedirs fast5s_native_1_23s_single     --alternate-fast5-basedirs fast5s_ivt_1_23s_single     --statistics-file-basename fast5s_native_1_23s_single     --store-p-value     --statistic-type ks --processes 50
