/*
 * Yanocomp compare to get mods
 */


process yanocompanalysis {

    publishDir "$params.outdir/yanofinal_16s", pattern: "*16s*", mode:'copy'
    publishDir "$params.outdir/yanofinal_23s", pattern: "*23s*", mode:'copy'

    tag "yanocomp compare to get mods"

    container '/home/gandalf/containers/yanocomp.sif'

    input:
    tuple val(rep), path (hdf5_native), path(hdf5_ivt)

    output:
    tuple path ("*.bed"), path ("*_sm_preds.json"), val (rep), emit: yano_out

    script:
    """
    yanocomp gmmtest -p 25 -n 1 --fdr-threshold 1 --min-ks 0 -c $hdf5_native -t $hdf5_ivt -o ${hdf5_native.simpleName}.bed -s ${hdf5_native.simpleName}_sm_preds.json
    """
}