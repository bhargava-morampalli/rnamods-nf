/*
 * Map fastqs to 16s rRNA reference using Minimap2 for native and IVT data
 */

process nanoget {

    publishDir "$params.outdir/nanostats_16s_native", pattern: "*native*16s*", mode:'copy'
    publishDir "$params.outdir/nanostats_16s_ivt", pattern: "*ivt*16s*", mode:'copy'
    publishDir "$params.outdir/nanostats_23s_native", pattern: "*native*23s*", mode:'copy'
    publishDir "$params.outdir/nanostats_23s_ivt", pattern: "*ivt*16s*", mode:'copy'
    
    tag "get read and alignment stats from bam"

    input:
    path bams

    output:
    path "*.feather", emit: nanostats

    script:

    """
    /home/bhargavam/nanoget/scripts/create_feather.py --threads 20 --bam $bams --output ${bams.simpleName}.feather
    """
}