/*
 * Yanocomp data prep step
 */


process yanocompprep {

    publishDir "$params.outdir/yano_${summarytext.simpleName}", mode:'copy'

    tag "prepare data for yanocomp"

    container '/home/bhargavam/Documents/containers/yanocomp.sif'

    input:
    path eventaligntext
    path summarytext

    output:
    path "*.hdf5", emit: yanohdf5
    val true, emit: yanoprepdone

    script:
    """
    yanocomp prep -p 50 -e $eventaligntext -s $summarytext -h ${summarytext.simpleName}.hdf5
    """
}