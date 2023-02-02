/*
 * Map fastqs to 16s rRNA reference using Minimap2 for native and IVT data
 */

process mapbam_16s {

    publishDir "$params.outdir/mappedbams_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/mappedbams_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    
    tag "map native and ivt fastqs to 16s rRNA reference using minimap2"

    container '/home/bhargavam/Documents/containers/minimap_samtools.sif'
    
    input:
    path reference_16s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: sortedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_16s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}