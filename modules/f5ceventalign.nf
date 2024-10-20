/*
 * Eventalign using f5c / nanopolish
 */

process f5ceventalign {

    publishDir "$params.outdir/eventalign_native_16s", pattern: "native*16s*", mode:'copy'
    publishDir "$params.outdir/eventalign_native_23s", pattern: "native*23s*", mode:'copy'
    publishDir "$params.outdir/eventalign_ivt_16s", pattern: "ivt*16s*", mode:'copy'
    publishDir "$params.outdir/eventalign_ivt_23s", pattern: "ivt*23s*", mode:'copy'

    tag "eventalign using f5c (nanopolish)"

    container '/home/bmorampa/containers/f5c_1.1--h0326b38_1.sif'

    input:
    tuple path(fastq), path(index), path(fai), path(gzi), path(readdb), path(bam), path(bai), val(rep), path(fast5s)
    path reference

    output:
    tuple path ("*_summary.txt"), path ("*_eventalign.txt"), val(rep), emit: tuple_ch1
    path "versions.yml", emit: versions

    script:
    """
    f5c eventalign --scale-events --signal-index --rna -r $fastq -b $bam -g $reference --summary ${bam.simpleName}_summary.txt --threads 20 > ${bam.simpleName}_eventalign.txt
    """
}