/*
 * Preprocess reads using nanodoc tool processing
 */

process nanodocpreprocess {

    publishDir "$params.outdir/${singlefast5s.simpleName}_ndocpreprocess", mode:'copy'

    tag "preprocess reads for nanodoc processing"

    input:
    path singlefast5s
    path reference

    output:
    path "*.pq"
    path "index.txt"
    val true, emit: ndocpre_done

    script:
    """
    source ~/nanoDoc/src/nanodoc/bin/activate
    python /home/bhargavam/nanoDoc/src/nanoDoc.py formatfile -i $singlefast5s -o ./ -r $reference -t 20
    deactivate
    """
}