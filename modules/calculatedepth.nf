/*
 * Calculate depths from 16s bams for native and IVT data
 */

process calculatedepth {

    publishDir "$params.outdir/depths_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/depths_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/depths_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/depths_23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "calculate depth info from sorted bam"

    container '/home/bhargavam/Documents/containers/samtools_1.16.1--h6899075_1.sif'

    input:
    path sortedbams

    output:
    path "*.txt", emit: depths

    script:

    """
    samtools depth -a -m 0 $sortedbams > ${sortedbams.simpleName}.txt
    """
}