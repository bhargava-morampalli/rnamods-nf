/*
 * Calculate depths from 16s bams for native and IVT data
 */

process extractmappedreads {

    publishDir "$params.outdir/mappedfastqs_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "extract mapped reads from sam files into fastq"

    container '/home/bhargavam/Documents/containers/samtools_1.16.1--h6899075_1.sif'

    input:
    path sams

    output:
    path "*.fastq", emit: mappedfastqs

    script:

    """
    samtools fastq -F4 $sams > ${sams.simpleName}.fastq
    """
}