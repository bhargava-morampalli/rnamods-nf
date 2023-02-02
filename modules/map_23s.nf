/*
 * Map fastqs to 23s rRNA reference using Minimap2 for native and IVT data
 */

process map_23s {

    publishDir "$params.outdir/sams_23s_native", pattern: "native*23s*", mode:'copy'
    publishDir "$params.outdir/sams_23s_ivt", pattern: "ivt*23s*", mode:'copy'
    
    tag "map native and ivt fastqs to 23s rRNA reference using minimap2"

    container '/home/bhargavam/Documents/containers/minimap2_2.24--h7132678_1.sif'

    input:
    path reference_23s
    path fastqs

    output:
    path "*_23s.sam", emit: sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_23s $fastqs > ${fastqs.simpleName}_23s.sam --secondary=no
    """
}