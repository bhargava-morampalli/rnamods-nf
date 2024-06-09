/*
 * Preprocess reads using nanodoc tool processing
 */

process nanodoccompare_23s {

    publishDir "$params.outdir/nanodocfinal_23s", mode:'copy'

    tag "preprocess reads for nanodoc processing"

    input:
    path reference
    path ndoc_ivt
    path ndoc_native
    val flag
    val flag

    output:
    path "*.txt"
    val true, emit: ndocfinal_done

    script:
    """
    source ~/nanoDoc/src/nanodoc/bin/activate
    python /home/bhargavam/nanoDoc/src/nanoDoc.py analysis -w /home/bhargavam/nanoDoc/weight5mer/ -p /home/bhargavam/nanoDoc/param20.txt -r $reference -rraw $ndoc_ivt -traw $ndoc_native -o ${ndoc_ivt.simpleName}.txt -s 1 -e 3370
    deactivate
    """
}