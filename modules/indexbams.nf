/*
 * Index sorted bams for native and IVT data using samtools index
 */

process indexbams {

    publishDir "$params.outdir/bams_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/bams_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/bams_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/bams_23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "index the sorted bam file"

    container '/home/bhargavam/Documents/containers/samtools_1.16.1--h6899075_1.sif'

    input:
    path sortedbams

    output:
    path "*.bam.bai", emit: sortedbamindex

    script:

    """
    samtools index $sortedbams
    """
}