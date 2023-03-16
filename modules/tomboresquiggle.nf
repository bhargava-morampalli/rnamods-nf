process tomboresquiggle {

    fair true
    publishDir "$params.outdir/${singlefast5s.simpleName}", mode:'copy'

    tag "resquiggle fast5s with tombo"

    container '/home/bhargavam/Documents/containers/tombo_new.sif'

    input:
    path singlefast5s
    path reference

    output:
    tuple path("*"), env(REP), optional: true, emit: resquiggledone_ch

    script:
    '''
    tombo resquiggle $singlefast5s $reference --rna --processes 50 --overwrite --num-most-common-errors 5
    REP=$(echo $path | cut -d '_' -f 2-)
    '''
}
