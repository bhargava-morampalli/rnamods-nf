/*
 * Convert multifast5s into single fast5s
 */

process multitosingle {

    publishDir "$params.outdir", mode:'copy'

    tag "convert multifast5s into single fast5s"

    container '/home/bhargavam/Documents/containers/ont-fast5-api_4.1.0--pyhdfd78af_0.sif'

    input:
    path multitosingleinputpath

    output:
    path "*", emit: singlefast5s_ch
    val true, emit:multitosingledone

    script:

    """
    mkdir -p ${multitosingleinputpath.simpleName}_single
    multi_to_single_fast5 --input_path $multitosingleinputpath --save_path ${multitosingleinputpath.simpleName}_single --recursive
    """
}