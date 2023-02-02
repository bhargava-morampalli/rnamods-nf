/*
 * Extract mapped read ids into text files from mapped fastq files - native and IVT data using seqkit
 */

process extract_readids_ivt {

    publishDir "$params.outdir/mappedfastas_ivt", pattern: "*.fasta", mode:'copy'
    publishDir "$params.outdir/mappedids_ivt", pattern: "*.txt", mode:'copy'

    tag "extract mapped read ids from fastqs to text file"

    container '/home/bhargavam/Documents/containers/seqkit_2.3.1--h9ee0642_0.sif'

    input:
    path mappedfastqs

    output:
    path "*.fasta", emit: mappedfastas
    path "*.txt", emit: mappedreadids
    
    script:

    """
    seqkit fq2fa $mappedfastqs -o ${mappedfastqs.simpleName}.fasta
    seqkit seq ${mappedfastqs.simpleName}.fasta -n -i > ${mappedfastqs.simpleName}.txt
    """
}