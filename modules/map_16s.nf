/*
 * Map fastqs to 16s rRNA reference using Minimap2 for native and IVT data
 */

process map_16s {

    publishDir "$params.outdir/sams_16s_native", pattern: "native*16s*", mode:'copy'
    publishDir "$params.outdir/sams_16s_ivt", pattern: "ivt*16s*", mode:'copy'
    
    tag "map native and ivt fastqs to 16s rRNA reference using minimap2"

    container '/home/bhargavam/Documents/containers/minimap2_2.24--h7132678_1.sif'

    input:
    path reference_16s
    path fastqs

    output:
    path "*_16s.sam", emit: sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_16s $fastqs > ${fastqs.simpleName}_16s.sam --secondary=no
    """
}