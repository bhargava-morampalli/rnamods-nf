/*
 * Dataprep for xpore
 */

process xporeanalysis {
    
    publishDir "$params.outdir/xporefinal_${yaml.simpleName}", mode:'copy'

    tag "xpore analysis step"

    container '/home/bhargavam/Documents/containers/xpore_2.1--pyh5e36f6f_0.sif'
    
    input:
    path yamlfile
    val flag
    val flag

    output:
    path "diffmod.table", emit: xporetable
    path "diffmod.log", emit: xporelog
    val true, emit: xporeanalysis

    script:
    """
    xpore diffmod --config $yamlfile
    """
}
