/*
 * index using f5c / nanopolish
 */

process f5cindex {

    publishDir "$params.outdir/mappedfastqs_16s_native", pattern: "native*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_native", pattern: "native*23s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_16s_ivt", pattern: "ivt*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_ivt", pattern: "ivt*23s*", mode:'copy'

    tag "nanopolish index using f5c"

    container '/home/bhargavam/Documents/containers/f5c_1.1--h0326b38_1.sif'

    input:
    path multifast5s
    path fastq

    output:
    path "*.fastq.index", emit: reads_index_ch
	path "*.fastq.index.fai", emit: reads_fai_ch
	path "*.fastq.index.gzi", emit: reads_gzi_ch
	path "*.fastq.index.readdb", emit: reads_readdb_ch
    val true, emit: f5cindexout

    script:
    """
    f5c index -t 24 -d $multifast5s $fastq
    """
}