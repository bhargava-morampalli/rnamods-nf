nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */


params.reference_16s = "$baseDir/references/k12_16s.fa"
params.reference_23s = "$baseDir/references/k12_23s.fa"

params.nativereads = "$baseDir/native/*.fastq"
params.ivtreads = "$baseDir/ivt/*.fastq"

params.nativefast5s = "/home/bmorampa/k12_native_fast5"
params.ivtfast5s = "/home/bmorampa/k12_ivt_fast5"


params.outdir = "results"


/*
 * pipeline processes
 */

include { map_16s as map_16s_native; map_16s as map_16s_ivt } from './modules/map_16s'

include { map_23s as map_23s_native; map_23s as map_23s_ivt } from './modules/map_23s'

include { mapstats_samtools as flagstat_16s_native; mapstats_samtools as flagstat_16s_ivt; mapstats_samtools as flagstat_23s_native; mapstats_samtools as flagstat_23s_ivt } from './modules/mapstats_samtools'

include { sortedbams as bam_16s_native; sortedbams as bam_16s_ivt; sortedbams as bam_23s_native; sortedbams as bam_23s_ivt } from './modules/sortedbams'

include { indexbams as bamindex_16s_native; indexbams as bamindex_16s_ivt; indexbams as bamindex_23s_native; indexbams as bamindex_23s_ivt } from './modules/indexbams'

include { calculatedepth as depths_16s_native; calculatedepth as depths_16s_ivt; calculatedepth as depths_23s_native; calculatedepth as depths_23s_ivt } from './modules/calculatedepth'

include { coverageplots as coverageplot_16s_native; coverageplots as coverageplot_16s_ivt; coverageplots as coverageplot_23s_native; coverageplots as coverageplot_23s_ivt } from './modules/coverageplots'

include { nanoget as nanostats_16s_native; nanoget as nanostats_16s_ivt; nanoget as nanostats_23s_native; nanoget as nanostats_23s_ivt } from './modules/nanoget'

include { extractmappedreads as extract_16s_native; extractmappedreads as extract_16s_ivt; extractmappedreads as extract_23s_native; extractmappedreads as extract_23s_ivt } from './modules/extractmappedreads'

include { extract_readids_native as nativereadidextract_16s; extract_readids_native as nativereadidextract_23s } from './modules/extract_readids_native'

include { extract_readids_ivt as ivtreadidextract_16s; extract_readids_ivt as ivtreadidextract_23s } from './modules/extract_readids_ivt'

include { extract_fast5s as extractfast5s_native_16s; extract_fast5s as extractfast5s_ivt_16s; extract_fast5s as extractfast5s_native_23s; extract_fast5s as extractfast5s_ivt_23s } from './modules/extract_fast5s'

include { mapbam_16s as mappedbams_16s_native; mapbam_16s as mappedbams_16s_ivt } from './modules/mapbam_16s'

include { mapbam_23s as mappedbams_23s_native; mapbam_23s as mappedbams_23s_ivt } from './modules/mapbam_23s'

include { indexmappedbams as mappedbamindex_16s_native; indexmappedbams as mappedbamindex_16s_ivt; indexmappedbams as mappedbamindex_23s_native; indexmappedbams as mappedbamindex_23s_ivt } from './modules/indexmappedbams'

include { multitosingle as multitosingle_n_16s;  multitosingle as multitosingle_n_23s } from './modules/multitosingle'

include { multitosingle as multitosingle_i_16s;  multitosingle as multitosingle_i_23s } from './modules/multitosingle'

include { tomboresquiggle as tomboresquiggle_n_16s; tomboresquiggle as tomboresquiggle_n_23s } from './modules/tomboresquiggle'

include { tomboresquiggle as tomboresquiggle_i_16s; tomboresquiggle as tomboresquiggle_i_23s } from './modules/tomboresquiggle'

include { tombocompare as tombocompare_16s; tombocompare as tombocompare_23s } from './modules/tombocompare'

include { tomboextract_16s as tomboextract_16s } from './modules/tomboextract_16s'

include { tomboextract_23s as tomboextract_23s } from './modules/tomboextract_23s'

include { f5cindex as f5cindex_n_16s; f5cindex as f5cindex_n_23s } from './modules/f5cindex'

include { f5cindex as f5cindex_i_16s; f5cindex as f5cindex_i_23s } from './modules/f5cindex'

include { f5ceventalign as f5ceventalign_n_16s; f5ceventalign as f5ceventalign_n_23s } from './modules/f5ceventalign'

include { f5ceventalign as f5ceventalign_i_16s; f5ceventalign as f5ceventalign_i_23s } from './modules/f5ceventalign'

include { yanocompprep as yanocompprep_n_16s; yanocompprep as yanocompprep_n_23s } from './modules/yanocompprep'

include { yanocompprep as yanocompprep_i_16s; yanocompprep as yanocompprep_i_23s } from './modules/yanocompprep'

include { yanocompanalysis as yanoanalysis_16s; yanocompanalysis as yanoanalysis_23s } from './modules/yanocompanalysis'

include { xporeprep as xporeprep_n_16s; xporeprep as xporeprep_n_23s } from './modules/xporeprep'

include { xporeprep as xporeprep_i_16s; xporeprep as xporeprep_i_23s } from './modules/xporeprep'


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
    
    tomboresquiggle_16s = tomboresquiggle_n_16s.out.resquiggledone_ch.mix(tomboresquiggle_i_16s.out.resquiggledone_ch).groupTuple(by: 1).view()
    tomboresquiggle_23s = tomboresquiggle_n_23s.out.resquiggledone_ch.mix(tomboresquiggle_i_23s.out.resquiggledone_ch).groupTuple(by: 1).view()

    tombocompare_16s(tomboresquiggle_16s)
    tombocompare_23s(tomboresquiggle_23s)

    tomboextract_16s (tombocompare_16s.out.tombostat_ch, tombocompare_16s.out.tombocomparedone)
    tomboextract_23s (tombocompare_23s.out.tombostat_ch, tombocompare_23s.out.tombocomparedone)

    f5c_index_n_16s = extractfast5s_native_16s.out.f5c_fast5.mix(extract_16s_native.out.f5c_fastq).groupTuple(by: 1).view()
    f5c_index_n_23s = extractfast5s_native_23s.out.f5c_fast5.mix(extract_23s_native.out.f5c_fastq).groupTuple(by: 1).view()

    f5c_index_i_16s = extractfast5s_ivt_16s.out.f5c_fast5.mix(extract_16s_ivt.out.f5c_fastq).groupTuple(by: 1).view()
    f5c_index_i_23s = extractfast5s_ivt_23s.out.f5c_fast5.mix(extract_23s_ivt.out.f5c_fastq).groupTuple(by: 1).view()

    f5cindex_n_16s(f5c_index_n_16s)
    f5cindex_n_23s(f5c_index_n_23s)

    f5cindex_i_16s(f5c_index_i_16s)
    f5cindex_i_23s(f5c_index_i_23s)

    f5cindex_n_16s.out.fastqindex.view { println "f5cindex output: $it" }
    mappedbamindex_16s_native.out.mappedbamindex.view { println "mappedbamindex output: $it" }

    f5c_event_n_16s = f5cindex_n_16s.out.fastqindex.mix(mappedbamindex_16s_native.out.mappedbamindex).groupTuple(by: 2).map { id, files, rep, fast5 -> tuple(files[1][0], files[1][1], files[1][2], files[1][3], files[1][4], files[0][0], files[0][1], rep, fast5) }.view()
    f5c_event_i_16s = f5cindex_i_16s.out.fastqindex.mix(mappedbamindex_16s_ivt.out.mappedbamindex).groupTuple(by: 2).map { id, files, rep, fast5 -> tuple(files[1][0], files[1][1], files[1][2], files[1][3], files[1][4], files[0][0], files[0][1], rep, fast5) }.view()
    f5c_event_n_23s = f5cindex_n_23s.out.fastqindex.mix(mappedbamindex_23s_native.out.mappedbamindex).groupTuple(by: 2).map { id, files, rep, fast5 -> tuple(files[1][0], files[1][1], files[1][2], files[1][3], files[1][4], files[0][0], files[0][1], rep, fast5) }.view()
    f5c_event_i_23s = f5cindex_i_23s.out.fastqindex.mix(mappedbamindex_23s_ivt.out.mappedbamindex).groupTuple(by: 2).map { id, files, rep, fast5 -> tuple(files[1][0], files[1][1], files[1][2], files[1][3], files[1][4], files[0][0], files[0][1], rep, fast5) }.view()

    f5ceventalign_n_16s(f5c_event_n_16s, reference_16s_ch)
    f5ceventalign_n_23s(f5c_event_n_23s, reference_23s_ch)
    
    f5ceventalign_i_16s(f5c_event_i_16s, reference_16s_ch)
    f5ceventalign_i_23s(f5c_event_i_23s, reference_23s_ch)

    yanocompprep_n_16s(f5ceventalign_n_16s.out.tuple_ch1)
    yanocompprep_n_23s(f5ceventalign_n_23s.out.tuple_ch1)
    
    yanocompprep_i_16s(f5ceventalign_i_16s.out.tuple_ch1)
    yanocompprep_i_23s(f5ceventalign_i_23s.out.tuple_ch1)

    yanocompprep_n_16s.out.yanohdf5.view()
    yanocompprep_i_16s.out.yanohdf5.view()

    yanocompprep_n_23s.out.yanohdf5.view()
    yanocompprep_i_23s.out.yanohdf5.view()

    yanoprep_16s = yanocompprep_n_16s.out.yanohdf5.join(yanocompprep_i_16s.out.yanohdf5, by:1).view { println "yanocompprep 16s output: $it" }
    yanoprep_23s = yanocompprep_n_23s.out.yanohdf5.join(yanocompprep_i_23s.out.yanohdf5, by:1).view { println "yanocompprep 23s output: $it" }

    yanoanalysis_16s(yanoprep_16s)
    yanoanalysis_23s(yanoprep_23s)

    yanoanalysis_16s.out.yano_out.view{ println "yanoanalysis 16s output: $it" }
    yanoanalysis_23s.out.yano_out.view{ println "yanoanalysis 23s output: $it" }

    xporeprep_n_16s(f5ceventalign_n_16s.out.tuple_ch1)
    xporeprep_n_23s(f5ceventalign_n_23s.out.tuple_ch1)

    xporeprep_i_16s(f5ceventalign_i_16s.out.tuple_ch1)
    xporeprep_i_23s(f5ceventalign_i_23s.out.tuple_ch1)

    xporeprep_n_16s.out.xporeout.join(xporeprep_i_16s.out.xporeout, by:1).view { println "xporeprep 16s output: $it" }
    xporeprep_n_23s.out.xporeout.join(xporeprep_i_23s.out.xporeout, by:1).view { println "xporeprep 23s output: $it" }

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
