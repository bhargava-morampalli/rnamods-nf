/*
 * Compare mods between native and ivt samples using tombo level sample compare
 */

process tombocompare {

    publishDir "$params.outdir/tombostat_16s", pattern: "*16s*", mode:'copy'
    publishDir "$params.outdir/tombostat_23s", pattern: "*23s*", mode:'copy'

    tag "tombo level sample compare"

    container '/home/bmorampa/containers/tombo_new.sif'

    input:
    tuple  path(singlefast5s), val(rep)

    output:
    path "*", emit: tombostat_ch
    val true, emit: tombocomparedone

    script:
    first = singlefast5s[0].baseName.toString()
    """
    type=\$(echo $first | cut -d '_' -f 2)
    echo \$type
    if [[ \$type == 'native' ]]; then 
        echo "native"
        tombo detect_modifications level_sample_compare \
        --fast5-basedirs ${singlefast5s[0]} \
        --alternate-fast5-basedirs ${singlefast5s[1]} \
        --statistics-file-basename ${singlefast5s[0].simpleName} \
        --store-p-value \
        --minimum-test-reads 1 \
        --statistic-type ks --processes 20

    else
        echo "ivt"
        tombo detect_modifications level_sample_compare \
        --fast5-basedirs ${singlefast5s[1]} \
        --alternate-fast5-basedirs ${singlefast5s[0]} \
        --statistics-file-basename ${singlefast5s[1].simpleName} \
        --store-p-value \
        --minimum-test-reads 1 \
        --statistic-type ks --processes 20
    fi
    """
}
