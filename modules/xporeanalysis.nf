/*
 * Dataprep for xpore
 */

process xporeanalysis {
    
    publishDir "$params.outdir/xporefinal_16s", pattern:"", mode:'copy'
    publishDir "$params.outdir/xporefinal_23s", pattern:"", mode:'copy'

    tag "xpore analysis step"

    container '/home/gandalf/containers/xpore_2.1--pyh5e36f6f_0.sif'
    
    input:
    tuple 

    output:
    path "diffmod*", emit: xpore_diffmod_outputs

    script:
    """
    xpore diffmod --config $yamlfile
    """
}
