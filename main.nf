nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */


params.reference_16s = "/home/bhargavam/Documents/rnamods-nf/references/k12_16s.fa"
params.reference_23s = "/home/bhargavam/Documents/rnamods-nf/references/k12_23s.fa"

params.nativereads = "$baseDir/native/*.fastq"
params.ivtreads = "$baseDir/ivt/*.fastq"

params.nativefast5s = "/home/bhargavam/Documents/data/k12_native_fast5_all"
params.ivtfast5s = "/home/bhargavam/Documents/data/k12_ivt_fast5_all"


params.outdir = "results"


/*
 * pipeline processes
 */

include { map_16s as map_16s_native; map_16s as map_16s_ivt } from '/home/bhargavam/Documents/nextflowmodules/map_16s'

include { map_23s as map_23s_native; map_23s as map_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/map_23s'

include { mapstats_samtools as flagstat_16s_native; mapstats_samtools as flagstat_16s_ivt; mapstats_samtools as flagstat_23s_native; mapstats_samtools as flagstat_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/mapstats_samtools'

include { sortedbams as bam_16s_native; sortedbams as bam_16s_ivt; sortedbams as bam_23s_native; sortedbams as bam_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/sortedbams'

include { indexbams as bamindex_16s_native; indexbams as bamindex_16s_ivt; indexbams as bamindex_23s_native; indexbams as bamindex_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/indexbams'

include { calculatedepth as depths_16s_native; calculatedepth as depths_16s_ivt; calculatedepth as depths_23s_native; calculatedepth as depths_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/calculatedepth'

include { coverageplots as coverageplot_16s_native; coverageplots as coverageplot_16s_ivt; coverageplots as coverageplot_23s_native; coverageplots as coverageplot_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/coverageplots'

include { nanoget as nanostats_16s_native; nanoget as nanostats_16s_ivt; nanoget as nanostats_23s_native; nanoget as nanostats_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/nanoget'

include { extractmappedreads as extract_16s_native; extractmappedreads as extract_16s_ivt; extractmappedreads as extract_23s_native; extractmappedreads as extract_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/extractmappedreads'

include { extract_readids_native as nativereadidextract_16s; extract_readids_native as nativereadidextract_23s } from '/home/bhargavam/Documents/nextflowmodules/extract_readids_native'

include { extract_readids_ivt as ivtreadidextract_16s; extract_readids_ivt as ivtreadidextract_23s } from '/home/bhargavam/Documents/nextflowmodules/extract_readids_ivt'

include { extract_fast5s as extractfast5s_native_16s; extract_fast5s as extractfast5s_ivt_16s; extract_fast5s as extractfast5s_native_23s; extract_fast5s as extractfast5s_ivt_23s } from '/home/bhargavam/Documents/nextflowmodules/extract_fast5s'

include { mapbam_16s as mappedbams_16s_native; mapbam_16s as mappedbams_16s_ivt } from '/home/bhargavam/Documents/nextflowmodules/mapbam_16s'

include { mapbam_23s as mappedbams_23s_native; mapbam_23s as mappedbams_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/mapbam_23s'

include { indexmappedbams as mappedbamindex_16s_native; indexmappedbams as mappedbamindex_16s_ivt; indexmappedbams as mappedbamindex_23s_native; indexmappedbams as mappedbamindex_23s_ivt } from '/home/bhargavam/Documents/nextflowmodules/indexmappedbams'

include { multitosingle as multitosingle_n_16s;  multitosingle as multitosingle_n_23s } from '/home/bhargavam/Documents/nextflowmodules/multitosingle'

include { multitosingle as multitosingle_i_16s;  multitosingle as multitosingle_i_23s } from '/home/bhargavam/Documents/nextflowmodules/multitosingle'

include { tomboresquiggle as tomboresquiggle_n_16s; tomboresquiggle as tomboresquiggle_n_23s } from '/home/bhargavam/Documents/nextflowmodules/tomboresquiggle'

include { tomboresquiggle as tomboresquiggle_i_16s; tomboresquiggle as tomboresquiggle_i_23s } from '/home/bhargavam/Documents/nextflowmodules/tomboresquiggle'


/*
 *  WORKFLOW 
 */



workflow {


    Channel
    .value(file(params.reference_16s))
    .set {reference_16s_ch}

    Channel
    .value(file(params.reference_23s))
    .set {reference_23s_ch}

    Channel
    .fromPath(params.nativereads)
    .set {nativefastqs_ch}

    Channel
    .fromPath(params.ivtreads)
    .set {ivtfastqs_ch}

    Channel
    .value(file(params.nativefast5s))
    .set {nativefast5s_ch}

    Channel
    .value(file(params.ivtfast5s))
    .set {ivtfast5s_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_native/*.fastq"))
    .set {mappedfastqs_16s_native}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_native/*.fastq"))
    .set {mappedfastqs_23s_native}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_ivt/*.fastq"))
    .set {mappedfastqs_16s_ivt}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_ivt/*.fastq"))
    .set {mappedfastqs_23s_ivt}


    map_16s_native (reference_16s_ch, nativefastqs_ch)
    map_16s_ivt (reference_16s_ch, ivtfastqs_ch)
    map_23s_native (reference_23s_ch, nativefastqs_ch)
    map_23s_ivt (reference_23s_ch, ivtfastqs_ch)

    flagstat_16s_native (map_16s_native.out.sams)
    flagstat_16s_ivt (map_16s_ivt.out.sams)
    flagstat_23s_native (map_23s_native.out.sams)
    flagstat_23s_ivt (map_23s_ivt.out.sams)

    bam_16s_native (map_16s_native.out.sams)
    bam_16s_ivt (map_16s_ivt.out.sams)
    bam_23s_native (map_23s_native.out.sams)
    bam_23s_ivt (map_23s_ivt.out.sams)

    bamindex_16s_native (bam_16s_native.out.sortedbams)
    bamindex_16s_ivt (bam_16s_ivt.out.sortedbams)
    bamindex_23s_native (bam_23s_native.out.sortedbams)
    bamindex_23s_ivt (bam_23s_ivt.out.sortedbams)

    depths_16s_native (bam_16s_native.out.sortedbams)
    depths_16s_ivt (bam_16s_ivt.out.sortedbams)
    depths_23s_native (bam_23s_native.out.sortedbams)
    depths_23s_ivt (bam_23s_ivt.out.sortedbams)

    coverageplot_16s_native (depths_16s_native.out.depths)
    coverageplot_16s_ivt (depths_16s_ivt.out.depths)
    coverageplot_23s_native (depths_23s_native.out.depths)
    coverageplot_23s_ivt (depths_23s_ivt.out.depths)

    nanostats_16s_native (bam_16s_native.out.sortedbams)
    nanostats_16s_ivt (bam_16s_ivt.out.sortedbams)
    nanostats_23s_native (bam_23s_native.out.sortedbams)
    nanostats_23s_ivt (bam_23s_ivt.out.sortedbams)

    extract_16s_native(map_16s_native.out.sams)
    extract_16s_ivt(map_16s_ivt.out.sams)
    extract_23s_native(map_23s_native.out.sams)
    extract_23s_ivt(map_23s_ivt.out.sams)

    nativereadidextract_16s(extract_16s_native.out.mappedfastqs)
    ivtreadidextract_16s(extract_16s_ivt.out.mappedfastqs)
    nativereadidextract_23s(extract_23s_native.out.mappedfastqs)
    ivtreadidextract_23s(extract_23s_ivt.out.mappedfastqs)

    extractfast5s_native_16s (nativefast5s_ch, nativereadidextract_16s.out.mappedreadids)
    extractfast5s_ivt_16s (ivtfast5s_ch, ivtreadidextract_16s.out.mappedreadids)
    extractfast5s_native_23s (nativefast5s_ch, nativereadidextract_23s.out.mappedreadids)
    extractfast5s_ivt_23s (ivtfast5s_ch, ivtreadidextract_23s.out.mappedreadids)

    mappedbams_16s_native (reference_16s_ch, extract_16s_native.out.mappedfastqs)
    mappedbams_16s_ivt (reference_16s_ch, extract_16s_ivt.out.mappedfastqs)
    mappedbams_23s_native (reference_23s_ch, extract_23s_native.out.mappedfastqs)
    mappedbams_23s_ivt (reference_23s_ch, extract_23s_ivt.out.mappedfastqs)

    mappedbamindex_16s_native (mappedbams_16s_native.out.sortedbams)
    mappedbamindex_16s_ivt (mappedbams_16s_ivt.out.sortedbams)
    mappedbamindex_23s_native (mappedbams_23s_native.out.sortedbams)
    mappedbamindex_23s_ivt (mappedbams_23s_ivt.out.sortedbams)


    multitosingle_n_16s (extractfast5s_native_16s.out.subsetfast5s)
    multitosingle_n_23s (extractfast5s_native_23s.out.subsetfast5s)
    multitosingle_i_16s (extractfast5s_ivt_16s.out.subsetfast5s)
    multitosingle_i_23s (extractfast5s_ivt_23s.out.subsetfast5s)

    tomboresquiggle_n_16s (multitosingle_n_16s.out.singlefast5s_ch, reference_16s_ch)
    tomboresquiggle_n_23s (multitosingle_n_23s.out.singlefast5s_ch, reference_23s_ch)
    tomboresquiggle_i_16s (multitosingle_i_16s.out.singlefast5s_ch, reference_16s_ch)
    tomboresquiggle_i_23s (multitosingle_i_23s.out.singlefast5s_ch, reference_23s_ch)

    tomboresquiggle_n_16s.out.resquiggledone
        .combine(tomboresquiggle_i_16s.out.resquiggledone)
        .view()
    
    left = Channel.of(['A', 1], ['B', 2], ['A', 3])
    right = Channel.of(['B', 'x'], ['B', 'y'], ['A', 'z'], ['A', 'w'])

    left
        .combine(right, by: 0)
        .view()

}


workflow.onComplete {

"""\
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    .stripIndent()

}