/*
 * Yanocomp data prep step
 */


process yanocompprep {

    publishDir "$params.outdir/yano_${summary.simpleName}", mode:'copy'

    tag "prepare data for yanocomp"

    container '/home/bmorampa/containers/yanocomp.sif'

    input:
    tuple path (summary), path (eventalign), val(rep)

    output:
    tuple path ("*.hdf5"), val (rep), emit: yanohdf5
    path "versions.yml", emit: versions

    script:
    """
    yanocomp prep -p 50 -e $eventalign -s $summary -h ${summary.simpleName}.hdf5
    """
}