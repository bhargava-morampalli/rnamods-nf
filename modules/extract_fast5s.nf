/*
 * Extract specific fast5s for all mapped reads for 16s and 23s - native and IVT data using fast5_subset
 */

process extract_fast5s {

    publishDir "$params.outdir", mode:'copy'
    
    tag "extract multifast5s for mapped reads"

    container '/home/bmorampa/containers/ont-fast5-api_4.1.0--pyhdfd78af_0.sif'

    input:
    path fast5inputpath
    path idtextfile

    output:
    path "fast5s*", emit: subsetfast5s
    val true, emit: extractdone_ch
    tuple path ("fast5s*"), env(REP), optional: true, emit: f5c_fast5
    path "versions.yml", emit: versions

    script:

    """
    mkdir -p fast5s_${idtextfile.simpleName}
    fast5_subset -i $fast5inputpath -s fast5s_${idtextfile.simpleName} -l *.txt -f ${idtextfile.simpleName}- --recursive
    REP=\$(echo fast5s_${idtextfile.simpleName} | cut -d '_' -f 3)
    """
}