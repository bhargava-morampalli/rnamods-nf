/*
 * Dataprep for xpore
 */


process xporeprep {

    publishDir "$params.outdir/xpore_${eventalign.simpleName}", mode:'copy'

    tag "prepare the data for xpore analysis"

    container '/home/bmorampa/containers/xpore_2.1--pyh5e36f6f_0.sif'

    input:
    tuple path(summary), path(eventalign), val(rep)

    output:
    tuple path("dataprep"), val(rep), emit: xporeout

    script:
    """
    xpore dataprep \
    --eventalign $eventalign \
    --out_dir dataprep \
    --n_processes 25
    """
}
