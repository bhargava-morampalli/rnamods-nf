/*
 * Yanocomp compare to get mods
 */


process yanocompanalysis {

    publishDir "$params.outdir/yanofinal_16s", pattern: "*16s*", mode:'copy'
    publishDir "$params.outdir/yanofinal_23s", pattern: "*23s*", mode:'copy'

    tag "yanocomp compare to get mods"

    container '/home/bhargavam/Documents/containers/yanocomp.sif'


    input:
    path hdf5_native
    path hdf5_ivt
    val flag
    val flag

    output:
    path "*.bed", emit: yanobed
    path "*_sm_preds.json", emit: yanojson
    val true, emit: yanocompanalysis

    script:
    """
    yanocomp gmmtest -p 40 --fdr-threshold 1 --min-ks 0 -c $hdf5_native -t $hdf5_ivt -o ${hdf5_native.simpleName}.bed -s ${hdf5_native.simpleName}_sm_preds.json
    """
}