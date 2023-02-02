/*
 * Compare mods between native and ivt samples using tombo level sample compare
 */

process tombocompare {

    publishDir "$params.outdir/tombostat_16s", pattern: "*16s*", mode:'copy'
    publishDir "$params.outdir/tombostat_23s", pattern: "*23s*", mode:'copy'

    tag "tombo level sample compare"

    container '/home/bhargavam/Documents/containers/tombo_new.sif'

    input:
    tuple path(native_singlefast5s), path (ivt_singlefast5s)
    val flag
    val flag

    output:
    path "*", emit: tombostat_ch
    val true, emit: tombocomparedone

    script:
    """
    tombo detect_modifications level_sample_compare \
    --fast5-basedirs $native_singlefast5s \
    --alternate-fast5-basedirs $ivt_singlefast5s \
    --statistics-file-basename ${native_singlefast5s.simpleName} \
    --store-p-value \
    --statistic-type ks --processes 50

    """
}