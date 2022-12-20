nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */


params.reference_16s = "/data/bhargava/reference_files/K12/k12_16s_88_extended.fa"
params.reference_23s = "/data/bhargava/reference_files/K12/k12_23s_78_extended.fa"

params.nativereads = "$baseDir/native/*.fastq"
params.ivtreads = "$baseDir/ivt/*.fastq"

params.nativefast5s = "/data/bhargava/k12_native_fast5/"
params.ivtfast5s = "/data/bhargava/k12_ivt_fast5/"


params.outdir = "results"


/*
 * pipeline processes
 */

/*
 * Get raw data from bam files for 16s rRNA - native and IVT data using nanoget module
 */

process nanostats_16s_native {

    publishDir "$params.outdir/nanostats_16s_native", mode:'copy'
    
    tag "get read and alignment stats from bam"

    input:
    path native16bams

    output:
    path "*.feather", emit: native16nanostats

    script:

    """
    /home/bhargavam/nanoget/scripts/create_feather.py --threads 20 --bam $native16bams --output ${native16bams.simpleName}.feather
    """
}

process nanostats_16s_ivt {

    publishDir "$params.outdir/nanostats_16s_ivt", mode:'copy'
    
    tag "get read and alignment stats from bam"

    input:
    path ivt16bams

    output:
    path "*.feather", emit: ivt16nanostats

    script:

    """
    /home/bhargavam/nanoget/scripts/create_feather.py --threads 20 --bam $ivt16bams --output ${ivt16bams.simpleName}.feather
    """
}

/*
 * Get raw data from bam files for 23s rRNA - native and IVT data using nanoget module
 */

process nanostats_23s_native {

    publishDir "$params.outdir/nanostats_23s_native", mode:'copy'
    
    tag "get read and alignment stats from bam"

    input:
    path native23bams

    output:
    path "*.feather", emit: native23nanostats

    script:

    """
    /home/bhargavam/nanoget/scripts/create_feather.py --threads 20 --bam $native23bams --output ${native23bams.simpleName}.feather
    """
}

process nanostats_23s_ivt {

    publishDir "$params.outdir/nanostats_23s_ivt", mode:'copy'
    
    tag "get read and alignment stats from bam"

    input:
    path ivt23bams

    output:
    path "*.feather", emit: ivt23nanostats
    
    script:

    """
    /home/bhargavam/nanoget/scripts/create_feather.py --threads 20 --bam $ivt23bams --output ${ivt23bams.simpleName}.feather
    """
}


include { extract_readids_native as nativereadidextract_16s; extract_readids_native as nativereadidextract_23s } from '/home/bhargavam/modules_nextflow/extract_readids_native'

include { extract_readids_ivt as ivtreadidextract_16s; extract_readids_ivt as ivtreadidextract_23s } from '/home/bhargavam/modules_nextflow/extract_readids_ivt'

include { extract_fast5s as extractfast5s_native_16s; extract_fast5s as extractfast5s_ivt_16s; extract_fast5s as extractfast5s_native_23s; extract_fast5s as extractfast5s_ivt_23s } from '/home/bhargavam/modules_nextflow/extract_fast5s'



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

    flagstat_16s_native (map_16s_native.out.native16sams)
    flagstat_16s_ivt (map_16s_ivt.out.ivt16sams)
    flagstat_23s_native (map_23s_native.out.native23sams)
    flagstat_23s_ivt (map_23s_ivt.out.ivt23sams)

    bam_16s_native (map_16s_native.out.native16sams)
    bam_16s_ivt (map_16s_ivt.out.ivt16sams)
    bam_23s_native (map_23s_native.out.native23sams)
    bam_23s_ivt (map_23s_ivt.out.ivt23sams)

    bamindex_16s_native (bam_16s_native.out.native16bams)
    bamindex_16s_ivt (bam_16s_ivt.out.ivt16bams)
    bamindex_23s_native (bam_23s_native.out.native23bams)
    bamindex_23s_ivt (bam_23s_ivt.out.ivt23bams)

    depths_16s_native (bam_16s_native.out.native16bams)
    depths_16s_ivt (bam_16s_ivt.out.ivt16bams)
    depths_23s_native (bam_23s_native.out.native23bams)
    depths_23s_ivt (bam_23s_ivt.out.ivt23bams)

    coverageplot_16s_native (depths_16s_native.out.native16depths)
    coverageplot_16s_ivt (depths_16s_ivt.out.ivt16depths)
    coverageplot_23s_native (depths_23s_native.out.native23depths)
    coverageplot_23s_ivt (depths_23s_ivt.out.ivt23depths)

    nanostats_16s_native (bam_16s_native.out.native16bams)
    nanostats_16s_ivt (bam_16s_ivt.out.ivt16bams)
    nanostats_23s_native (bam_23s_native.out.native23bams)
    nanostats_23s_ivt (bam_23s_ivt.out.ivt23bams)

    extract_16s_native(map_16s_native.out.native16sams)
    extract_16s_ivt(map_16s_ivt.out.ivt16sams)
    extract_23s_native(map_23s_native.out.native23sams)
    extract_23s_ivt(map_23s_ivt.out.ivt23sams)

    nativereadidextract_16s(extract_16s_native.out.native16mappedfastqs)
    ivtreadidextract_16s(extract_16s_ivt.out.ivt16mappedfastqs)
    nativereadidextract_23s(extract_23s_native.out.native23mappedfastqs)
    ivtreadidextract_23s(extract_23s_ivt.out.ivt23mappedfastqs)

    extractfast5s_native_16s (nativefast5s_ch, nativereadidextract_16s.out.mappedreadids)
    extractfast5s_ivt_16s (ivtfast5s_ch, ivtreadidextract_16s.out.mappedreadids)
    extractfast5s_native_23s (nativefast5s_ch, nativereadidextract_23s.out.mappedreadids)
    extractfast5s_ivt_23s (ivtfast5s_ch, ivtreadidextract_23s.out.mappedreadids)

    mappedbams_16s_native (reference_16s_ch, extract_16s_native.out.native16mappedfastqs)
    mappedbams_16s_ivt (reference_16s_ch, extract_16s_ivt.out.ivt16mappedfastqs)
    mappedbams_23s_native (reference_23s_ch, extract_23s_native.out.native23mappedfastqs)
    mappedbams_23s_ivt (reference_23s_ch, extract_23s_ivt.out.ivt23mappedfastqs)

    mappedbamindex_16s_native (mappedbams_16s_native.out.native16mappedbams)
    mappedbamindex_16s_ivt (mappedbams_16s_ivt.out.ivt16mappedbams)
    mappedbamindex_23s_native (mappedbams_23s_native.out.native23mappedbams)
    mappedbamindex_23s_ivt (mappedbams_23s_ivt.out.ivt23mappedbams)



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