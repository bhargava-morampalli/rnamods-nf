/*
 * Map fastqs to 23s rRNA reference using Minimap2 for native and IVT data
 */

process mapbam_23s {

    publishDir "$params.outdir/mappedbams_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/mappedbams_23s_ivt", pattern: "*ivt*23s*", mode:'copy'
    
    tag "map native and ivt fastqs to 23s rRNA reference using minimap2"

    container '/home/bhargavam/Documents/containers/minimap_samtools.sif'

    input:
    path reference_23s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: sortedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_23s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}