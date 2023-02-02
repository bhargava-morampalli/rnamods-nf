/*
 * Extract raw data from tombo stats files
 */

process tomboextract_16s {

    publishDir "$params.outdir/tombofinal_16s", mode:'copy'

    tag "extracting p values from tombo stat files for 16s"

    container '/home/bhargavam/Documents/containers/tombo_new.sif'

    input:
    path statfile
    val flag

    output:
    path "*.csv", emit: tomboextract_ch
    val true, emit: tomboextractdone

    script:
    """
    #! /usr/bin/env python
    from tombo import tombo_helper, tombo_stats, resquiggle
    import pandas as pd

    sample_level_stats = tombo_stats.LevelStats("$statfile")
    reg_level_stats = sample_level_stats.get_reg_stats('16s_88_rrsE', '+', 1, 1813)
    pd.DataFrame(reg_level_stats).to_csv("${statfile.simpleName}.csv")

    """
}
