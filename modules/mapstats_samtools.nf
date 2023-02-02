/*
 * get mapping statistics for all sam files using samtools flagstat
 */

process mapstats_samtools {

    publishDir "$params.outdir/mapstats_flagstat", mode:'copy'
    
    tag "mapping stats with samtools flagstat"

    container '/home/bhargavam/Documents/containers/samtools_1.16.1--h6899075_1.sif'

    input:
    path sams

    output:
    path "*_flagstat.txt", emit: flagstatfile

    script:

    """
    samtools flagstat $sams | tee ${sams.simpleName}_flagstat.txt
    """
}