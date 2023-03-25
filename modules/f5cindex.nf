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
    tuple path(fast5_fastq), val(rep)

    output:
    tuple path("*.fastq*")

    script:
    """
    type=\$(echo $first | cut -d '.' -f 2)
    echo \$type
    if [[ \$type == 'fastq' ]]; then
    f5c index -t 24 -d ${fast5_fastq[1]} ${fast5_fastq[0]}
    else
    f5c index -t 24 -d ${fast5_fastq[0]} ${fast5_fastq[1]}
    fi
    """
}