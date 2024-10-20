/*
 * Index mapped sorted bams for native and IVT data using samtools index
 */

process indexmappedbams {

    publishDir "$params.outdir/mappedbams_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/mappedbams_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/mappedbams_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/mappedbams_23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "index the mapped sorted bam file"

    container '/home/bmorampa/containers/samtools_1.16.1--h6899075_1.sif'

    input:
    path mappedsortedbams

    output:
    tuple val(mappedsortedbams), path ("*.bam*", includeInputs:true), env(REP), emit: mappedbamindex
    path "versions.yml", emit: versions

    script:

    """
    samtools index $mappedsortedbams
    REP=\$(echo ${mappedsortedbams.simpleName} | cut -d '_' -f 2)
    """
}