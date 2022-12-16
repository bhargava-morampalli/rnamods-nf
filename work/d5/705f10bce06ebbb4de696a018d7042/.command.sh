#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
file = pd.read_csv("native_1_23s_sorted.txt", sep='	', header=None)
file.columns = ["chrom", "position", "depth"]
sns.set_theme(style="darkgrid")
x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
plt.title("native - 23s rRNA")
plt.xlabel("position on 23s rRNA")
plt.ylabel("Coverage")
plt.savefig("native_1_23s_sorted" + ".pdf")
