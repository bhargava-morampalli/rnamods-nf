/*
 * Get coverage plots from depth files for native and IVT data
 */

process coverageplots {

    publishDir "$params.outdir/coverageplots16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/coverageplots16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/coverageplots23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/coverageplots23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "coverage plot from depth file"

    input:
    path depths

    output:
    path "*.pdf", emit: coverageplot

    script:

    """
    #! /home/bhargavam/mambaforge/bin/python
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    file = pd.read_csv("$depths", sep='\t', header=None)
    file.columns = ["chrom", "position", "depth"]
    sns.set_theme(style="darkgrid")
    x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
    plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
    plt.title("native - 16s rRNA")
    plt.xlabel("position on 16s rRNA")
    plt.ylabel("Coverage")
    plt.savefig("${depths.simpleName}" + ".pdf")
    """
}