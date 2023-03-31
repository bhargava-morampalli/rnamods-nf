/*
 * Eventalign using f5c / nanopolish
 */

process f5ceventalign {

    publishDir "$params.outdir/eventalign_native_16s", pattern: "native*16s*", mode:'copy'
    publishDir "$params.outdir/eventalign_native_23s", pattern: "native*23s*", mode:'copy'
    publishDir "$params.outdir/eventalign_ivt_16s", pattern: "ivt*16s*", mode:'copy'
    publishDir "$params.outdir/eventalign_ivt_23s", pattern: "ivt*23s*", mode:'copy'

    tag "eventalign using f5c (nanopolish)"

    container '/home/bhargavam/Documents/containers/f5c_1.1--h0326b38_1.sif'

    input:
    tuple path(bam_fastq), val(rep)
    path reference

    output:
    path "*_summary.txt", emit: f5csummary_ch
    path "*_eventalign.txt", emit: f5ceventalign_ch
    val true, emit: f5ceventaligndone_ch

    script:
    firsttuple = bam_fastq[0]
    bam = firsttuple[0]
    fastq = firsttuple[1]
    """
    f5c eventalign --scale-events --signal-index --print-read-names --rna -r $fastq -b $bam -g $reference --summary ${bam_fastq[1].simpleName}_summary.txt --threads 40 > ${bam_fastq[1].simpleName}_eventalign.txt
    """
}