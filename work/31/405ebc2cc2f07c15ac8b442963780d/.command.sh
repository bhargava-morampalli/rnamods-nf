#!/usr/bin/env python
from tombo import tombo_helper, tombo_stats, resquiggle
import pandas as pd

sample_level_stats = tombo_stats.LevelStats("fast5s_native_3_23s_single.tombo.stats")
reg_level_stats = sample_level_stats.get_reg_stats('23s1_extended', '+', 1, 3163)
pd.DataFrame(reg_level_stats).to_csv("fast5s_native_3_23s_single.csv")
