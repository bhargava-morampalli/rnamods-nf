/*
 * Dataprep for xpore
 */


process xporeprep {

    publishDir "$params.outdir/xpore_${eventalign.simpleName}", pattern: 'dataprep/*', mode:'copy'

    tag "prepare the data for xpore analysis"

    container '/home/gandalf/containers/xpore_2.1--pyh5e36f6f_0.sif'

    input:
    tuple path(summary), path(eventalign), val(rep)

    output:
    path "dataprep/data.log", emit: xporelog
    path "dataprep/data.index", emit: xporeindex
    path "dataprep/data.json", emit: xporejson
    path "dataprep/data.readcount", emit: xporereadcount
    path "dataprep/eventalign.index", emit: eventalignindex
    val rep, emit: xporeprepdone

    script:
    """
    xpore dataprep \
    --eventalign $eventalign \
    --out_dir dataprep \
    --n_processes 25
    """
}
