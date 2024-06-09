process xporeroc16 {

    publishDir "$params.outdir/xporeroc23files", mode:'copy'

    input:
    path resultstable

    output:
    path "*.csv", emit: xporeroc16csv
    val true, emit: xporeroc16done

    script:
    """
    python3 /home/bhargavam/scripts/xpore_roc_23s.py $resultstable 1
    """
}