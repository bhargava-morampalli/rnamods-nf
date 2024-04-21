/*
 * index using f5c / nanopolish
 */

process f5cindex {

    publishDir "$params.outdir/mappedfastqs_16s_native", pattern: "native*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_native", pattern: "native*23s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_16s_ivt", pattern: "ivt*16s*", mode:'copy'
    publishDir "$params.outdir/mappedfastqs_23s_ivt", pattern: "ivt*23s*", mode:'copy'

    tag "nanopolish index using f5c"

    container '/home/gandalf/containers/f5c_1.1--h0326b38_1.sif'

    input:
    tuple path(fast5_fastq), val(rep)

    output:
    tuple val(second), path("*.fastq*"), val(rep), emit: fastqindex

    script:
    if (!fast5_fastq[0].baseName.contains("fast5")) {
        second = fast5_fastq[0]
    } else {
        second = fast5_fastq[1]
    }
    first = fast5_fastq[0].baseName.toString()
    """
    type=\$(echo $first | cut -d '_' -f 1)
    echo \$type
    if [[ \$type == 'fast5' ]]; then
    f5c index -t 24 -d ${fast5_fastq[0]} ${fast5_fastq[1]}
    else
    f5c index -t 24 -d ${fast5_fastq[1]} ${fast5_fastq[0]}
    fi
    """
}