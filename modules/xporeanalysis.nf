/*
 * Dataprep for xpore
 */

process xporeanalysis {
    
    publishDir "$params.outdir/xporefinal_16s", pattern:"", mode:'copy'
    publishDir "$params.outdir/xporefinal_23s", pattern:"", mode:'copy'

    tag "xpore analysis step"

    container '/home/bmorampa/containers/xpore_2.1--pyh5e36f6f_0.sif'
    
    input:
    tuple val(rep), path(xpore_native), path(xpore_ivt)

    output:
    path "diffmod*", emit: xpore_diffmod_outputs

    script:
    diffmod_config = "--config diffmod_config.yml"
    """
    create_yml.py diffmod_config.yml $xpore_native $xpore_ivt
    xpore diffmod $diffmod_config --n_processes $task.cpus
    """
}
