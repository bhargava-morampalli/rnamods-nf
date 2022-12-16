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


include { multitosingle as multitosingle_n_1_16s; multitosingle as multitosingle_n_2_16s; multitosingle as multitosingle_n_3_16s; multitosingle as multitosingle_n_1_23s; multitosingle as multitosingle_n_2_23s; multitosingle as multitosingle_n_3_23s } from '/home/bhargavam/modules_nextflow/multitosingle'

include { multitosingle as multitosingle_i_1_16s; multitosingle as multitosingle_i_2_16s; multitosingle as multitosingle_i_3_16s; multitosingle as multitosingle_i_1_23s; multitosingle as multitosingle_i_2_23s; multitosingle as multitosingle_i_3_23s } from '/home/bhargavam/modules_nextflow/multitosingle'

include { tomboresquiggle as tomboresquiggle_n_1_16s; tomboresquiggle as tomboresquiggle_n_2_16s; tomboresquiggle as tomboresquiggle_n_3_16s; tomboresquiggle as tomboresquiggle_n_1_23s; tomboresquiggle as tomboresquiggle_n_2_23s; tomboresquiggle as tomboresquiggle_n_3_23s } from '/home/bhargavam/modules_nextflow/tomboresquiggle'

include { tomboresquiggle as tomboresquiggle_i_1_16s; tomboresquiggle as tomboresquiggle_i_2_16s; tomboresquiggle as tomboresquiggle_i_3_16s; tomboresquiggle as tomboresquiggle_i_1_23s; tomboresquiggle as tomboresquiggle_i_2_23s; tomboresquiggle as tomboresquiggle_i_3_23s } from '/home/bhargavam/modules_nextflow/tomboresquiggle'

include { f5cindex as f5cindex1_n_16s; f5cindex as f5cindex1_n_23s ; f5cindex as f5cindex1_i_16s; f5cindex as f5cindex1_i_23s } from '/home/bhargavam/modules_nextflow/f5cindex'

include { f5cindex as f5cindex2_n_16s; f5cindex as f5cindex2_n_23s ; f5cindex as f5cindex2_i_16s; f5cindex as f5cindex2_i_23s} from '/home/bhargavam/modules_nextflow/f5cindex'

include { f5cindex as f5cindex3_n_16s; f5cindex as f5cindex3_n_23s ; f5cindex as f5cindex3_i_16s; f5cindex as f5cindex3_i_23s} from '/home/bhargavam/modules_nextflow/f5cindex'

include { f5ceventalign as f5ceventalign1_n_16s; f5ceventalign as f5ceventalign1_n_23s ; f5ceventalign as f5ceventalign1_i_16s; f5ceventalign as f5ceventalign1_i_23s} from '/home/bhargavam/modules_nextflow/f5ceventalign'

include { f5ceventalign as f5ceventalign2_n_16s; f5ceventalign as f5ceventalign2_n_23s ; f5ceventalign as f5ceventalign2_i_16s; f5ceventalign as f5ceventalign2_i_23s} from '/home/bhargavam/modules_nextflow/f5ceventalign'

include { f5ceventalign as f5ceventalign3_n_16s; f5ceventalign as f5ceventalign3_n_23s ; f5ceventalign as f5ceventalign3_i_16s; f5ceventalign as f5ceventalign3_i_23s} from '/home/bhargavam/modules_nextflow/f5ceventalign'

include { tombocompare as tombocompare_1_16s; tombocompare as tombocompare_2_16s; tombocompare as tombocompare_3_16s; tombocompare as tombocompare_1_23s; tombocompare as tombocompare_2_23s; tombocompare as tombocompare_3_23s } from '/home/bhargavam/modules_nextflow/tombocompare'

include { tomboextract_16s as tomboextract_16s_1 ; tomboextract_16s as tomboextract_16s_2 ; tomboextract_16s as tomboextract_16s_3 } from '/home/bhargavam/modules_nextflow/tomboextract_16s'

include { tomboextract_23s as tomboextract_23s_1 ; tomboextract_23s as tomboextract_23s_2 ; tomboextract_23s as tomboextract_23s_3 } from '/home/bhargavam/modules_nextflow/tomboextract_23s'


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
    .fromPath("$baseDir/results/fast5s_native_1_16s")
    .set {native_1_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_native_2_16s")
    .set {native_2_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_native_3_16s")
    .set {native_3_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_native_1_23s")
    .set {native_1_23s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_native_2_23s")
    .set {native_2_23s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_native_3_23s")
    .set {native_3_23s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_1_16s")
    .set {ivt_1_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_2_16s")
    .set {ivt_2_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_3_16s")
    .set {ivt_3_16s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_1_23s")
    .set {ivt_1_23s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_2_23s")
    .set {ivt_2_23s_fast5s_ch}

    Channel
    .fromPath("$baseDir/results/fast5s_ivt_3_23s")
    .set {ivt_3_23s_fast5s_ch}




    Channel
    .value(file("$baseDir/results/fast5s_native_1_16s_single"))
    .set {native_1_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_native_2_16s_single"))
    .set {native_2_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_native_3_16s_single"))
    .set {native_3_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_native_1_23s_single"))
    .set {native_1_23s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_native_2_23s_single"))
    .set {native_2_23s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_native_3_23s_single"))
    .set {native_3_23s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_1_16s_single"))
    .set {ivt_1_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_2_16s_single"))
    .set {ivt_2_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_3_16s_single"))
    .set {ivt_3_16s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_1_23s_single"))
    .set {ivt_1_23s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_2_23s_single"))
    .set {ivt_2_23s_single_fast5s_ch}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_3_23s_single"))
    .set {ivt_3_23s_single_fast5s_ch}





    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq"))
    .set {native_16s_1_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq"))
    .set {native_16s_2_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq"))
    .set {native_16s_3_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq"))
    .set {ivt_16s_1_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq"))
    .set {ivt_16s_2_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq"))
    .set {ivt_16s_3_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq"))
    .set {native_23s_1_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq"))
    .set {native_23s_2_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq"))
    .set {native_23s_3_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq"))
    .set {ivt_23s_1_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq"))
    .set {ivt_23s_2_fastqs_ch}

    Channel
    .value(file("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq"))
    .set {ivt_23s_3_fastqs_ch}


fastqinput1 = [
        ("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_native/native_1_16s.fastq.index.readdb")
]

fastqinput2 = [
        ("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_native/native_2_16s.fastq.index.readdb")
]

fastqinput3 = [
        ("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_native/native_3_16s.fastq.index.readdb")
]

fastqinput4 = [
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_1_16s.fastq.index.readdb")
]

fastqinput5 = [
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_2_16s.fastq.index.readdb")
]

fastqinput6 = [
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq.index"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_16s_ivt/ivt_3_16s.fastq.index.readdb")
]

fastqinput7 = [
        ("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_native/native_1_23s.fastq.index.readdb")
]

fastqinput8 = [
        ("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_native/native_2_23s.fastq.index.readdb")
]

fastqinput9 = [
        ("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_native/native_3_23s.fastq.index.readdb")
]

fastqinput10 = [
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_1_23s.fastq.index.readdb")
]

fastqinput11 = [
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_2_23s.fastq.index.readdb")
]

fastqinput12 = [
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq.index"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq.index.fai"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq.index.gzi"),
        ("$baseDir/results/mappedfastqs_23s_ivt/ivt_3_23s.fastq.index.readdb")
]



baminput1 = [
        file("$baseDir/results/mappedbams_16s_native/native_1_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_native/native_1_16s_sorted.bam.bai")
]

baminput2 = [
        file("$baseDir/results/mappedbams_16s_native/native_2_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_native/native_2_16s_sorted.bam.bai")
]

baminput3 = [
        file("$baseDir/results/mappedbams_16s_native/native_3_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_native/native_3_16s_sorted.bam.bai")
]

baminput4 = [
        file("$baseDir/results/mappedbams_16s_ivt/ivt_1_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_ivt/ivt_1_16s_sorted.bam.bai")
]

baminput5 = [
        file("$baseDir/results/mappedbams_16s_ivt/ivt_2_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_ivt/ivt_2_16s_sorted.bam.bai")
]

baminput6 = [
        file("$baseDir/results/mappedbams_16s_ivt/ivt_3_16s_sorted.bam"),
        file("$baseDir/results/mappedbams_16s_ivt/ivt_3_16s_sorted.bam.bai")
]

baminput7 = [
        file("$baseDir/results/mappedbams_23s_native/native_1_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_native/native_1_23s_sorted.bam.bai")
]

baminput8 = [
        file("$baseDir/results/mappedbams_23s_native/native_2_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_native/native_2_23s_sorted.bam.bai")
]

baminput9 = [
        file("$baseDir/results/mappedbams_23s_native/native_3_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_native/native_3_23s_sorted.bam.bai")
]

baminput10 = [
        file("$baseDir/results/mappedbams_23s_ivt/ivt_1_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_ivt/ivt_1_23s_sorted.bam.bai")
]

baminput11 = [
        file("$baseDir/results/mappedbams_23s_ivt/ivt_2_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_ivt/ivt_2_23s_sorted.bam.bai")
]

baminput12 = [
        file("$baseDir/results/mappedbams_23s_ivt/ivt_3_23s_sorted.bam"),
        file("$baseDir/results/mappedbams_23s_ivt/ivt_3_23s_sorted.bam.bai")
]





    multitosingle_n_1_16s (native_1_16s_fast5s_ch)
    multitosingle_n_2_16s (native_2_16s_fast5s_ch)
    multitosingle_n_3_16s (native_3_16s_fast5s_ch)

    multitosingle_n_1_23s (native_1_23s_fast5s_ch)
    multitosingle_n_2_23s (native_2_23s_fast5s_ch)
    multitosingle_n_3_23s (native_3_23s_fast5s_ch)

    multitosingle_i_1_16s (ivt_1_16s_fast5s_ch)
    multitosingle_i_2_16s (ivt_2_16s_fast5s_ch)
    multitosingle_i_3_16s (ivt_3_16s_fast5s_ch)

    multitosingle_i_1_23s (ivt_1_23s_fast5s_ch)
    multitosingle_i_2_23s (ivt_2_23s_fast5s_ch)
    multitosingle_i_3_23s (ivt_3_23s_fast5s_ch)

    tomboresquiggle_n_1_16s (native_1_16s_single_fast5s_ch, reference_16s_ch, multitosingle_n_1_16s.out.multitosingledone)
    tomboresquiggle_n_2_16s (native_2_16s_single_fast5s_ch, reference_16s_ch, multitosingle_n_2_16s.out.multitosingledone)
    tomboresquiggle_n_3_16s (native_3_16s_single_fast5s_ch, reference_16s_ch, multitosingle_n_3_16s.out.multitosingledone)

    tomboresquiggle_n_1_23s (native_1_23s_single_fast5s_ch, reference_23s_ch, multitosingle_n_1_23s.out.multitosingledone)
    tomboresquiggle_n_2_23s (native_2_23s_single_fast5s_ch, reference_23s_ch, multitosingle_n_2_23s.out.multitosingledone)
    tomboresquiggle_n_3_23s (native_3_23s_single_fast5s_ch, reference_23s_ch, multitosingle_n_3_23s.out.multitosingledone)

    tomboresquiggle_i_1_16s (ivt_1_16s_single_fast5s_ch, reference_16s_ch, multitosingle_i_1_16s.out.multitosingledone)
    tomboresquiggle_i_2_16s (ivt_2_16s_single_fast5s_ch, reference_16s_ch, multitosingle_i_2_16s.out.multitosingledone)
    tomboresquiggle_i_3_16s (ivt_3_16s_single_fast5s_ch, reference_16s_ch, multitosingle_i_3_16s.out.multitosingledone)

    tomboresquiggle_i_1_23s (ivt_1_23s_single_fast5s_ch, reference_23s_ch, multitosingle_i_1_23s.out.multitosingledone)
    tomboresquiggle_i_2_23s (ivt_2_23s_single_fast5s_ch, reference_23s_ch, multitosingle_i_2_23s.out.multitosingledone)
    tomboresquiggle_i_3_23s (ivt_3_23s_single_fast5s_ch, reference_23s_ch, multitosingle_i_3_23s.out.multitosingledone)

    tombocompare_1_16s (native_1_16s_single_fast5s_ch, ivt_1_16s_single_fast5s_ch, tomboresquiggle_n_1_16s.out.resquiggledone_ch, tomboresquiggle_i_1_16s.out.resquiggledone_ch)
    tombocompare_2_16s (native_2_16s_single_fast5s_ch, ivt_2_16s_single_fast5s_ch, tomboresquiggle_n_2_16s.out.resquiggledone_ch, tomboresquiggle_i_2_16s.out.resquiggledone_ch)
    tombocompare_3_16s (native_3_16s_single_fast5s_ch, ivt_3_16s_single_fast5s_ch, tomboresquiggle_n_3_16s.out.resquiggledone_ch, tomboresquiggle_i_3_16s.out.resquiggledone_ch)

    tombocompare_1_23s (native_1_23s_single_fast5s_ch, ivt_1_23s_single_fast5s_ch, tomboresquiggle_n_1_23s.out.resquiggledone_ch, tomboresquiggle_i_1_23s.out.resquiggledone_ch)
    tombocompare_2_23s (native_2_23s_single_fast5s_ch, ivt_2_23s_single_fast5s_ch, tomboresquiggle_n_2_23s.out.resquiggledone_ch, tomboresquiggle_i_2_23s.out.resquiggledone_ch)
    tombocompare_3_23s (native_3_23s_single_fast5s_ch, ivt_3_23s_single_fast5s_ch, tomboresquiggle_n_3_23s.out.resquiggledone_ch, tomboresquiggle_i_3_23s.out.resquiggledone_ch)


    f5cindex1_n_16s (native_1_16s_fast5s_ch, native_16s_1_fastqs_ch)
    f5cindex2_n_16s (native_2_16s_fast5s_ch, native_16s_2_fastqs_ch)
    f5cindex3_n_16s (native_3_16s_fast5s_ch, native_16s_3_fastqs_ch)

    f5cindex1_i_16s (ivt_1_16s_fast5s_ch, ivt_16s_1_fastqs_ch)
    f5cindex2_i_16s (ivt_2_16s_fast5s_ch, ivt_16s_2_fastqs_ch)
    f5cindex3_i_16s (ivt_3_16s_fast5s_ch, ivt_16s_3_fastqs_ch)

    f5cindex1_n_23s (native_1_23s_fast5s_ch, native_23s_1_fastqs_ch)
    f5cindex2_n_23s (native_2_23s_fast5s_ch, native_23s_2_fastqs_ch)
    f5cindex3_n_23s (native_3_23s_fast5s_ch, native_23s_3_fastqs_ch)

    f5cindex1_i_23s (ivt_1_23s_fast5s_ch, ivt_23s_1_fastqs_ch)
    f5cindex2_i_23s (ivt_2_23s_fast5s_ch, ivt_23s_2_fastqs_ch)
    f5cindex3_i_23s (ivt_3_23s_fast5s_ch, ivt_23s_3_fastqs_ch)

    f5ceventalign1_n_16s (fastqinput1, baminput1, reference_16s_ch, f5cindex1_n_16s.out.f5cindexout)
    f5ceventalign2_n_16s (fastqinput2, baminput2, reference_16s_ch, f5cindex2_n_16s.out.f5cindexout)
    f5ceventalign3_n_16s (fastqinput3, baminput3, reference_16s_ch, f5cindex3_n_16s.out.f5cindexout)

    f5ceventalign1_i_16s (fastqinput4, baminput4, reference_16s_ch, f5cindex1_i_16s.out.f5cindexout)
    f5ceventalign2_i_16s (fastqinput5, baminput5, reference_16s_ch, f5cindex2_i_16s.out.f5cindexout)
    f5ceventalign3_i_16s (fastqinput6, baminput6, reference_16s_ch, f5cindex3_i_16s.out.f5cindexout)

    f5ceventalign1_n_23s (fastqinput7, baminput7, reference_23s_ch, f5cindex1_n_23s.out.f5cindexout)
    f5ceventalign2_n_23s (fastqinput8, baminput8, reference_23s_ch, f5cindex2_n_23s.out.f5cindexout)
    f5ceventalign3_n_23s (fastqinput9, baminput9, reference_23s_ch, f5cindex3_n_23s.out.f5cindexout)

    f5ceventalign1_i_23s (fastqinput10, baminput10, reference_23s_ch, f5cindex1_i_23s.out.f5cindexout)
    f5ceventalign2_i_23s (fastqinput11, baminput11, reference_23s_ch, f5cindex2_i_23s.out.f5cindexout)
    f5ceventalign3_i_23s (fastqinput12, baminput12, reference_23s_ch, f5cindex3_i_23s.out.f5cindexout)

    tomboextract_16s_1 (tombocompare_1_16s.out.tombostat_ch, tombocompare_1_16s.out.tombocomparedone)
    tomboextract_16s_2 (tombocompare_2_16s.out.tombostat_ch, tombocompare_2_16s.out.tombocomparedone)
    tomboextract_16s_3 (tombocompare_3_16s.out.tombostat_ch, tombocompare_3_16s.out.tombocomparedone)

    tomboextract_23s_1 (tombocompare_1_23s.out.tombostat_ch, tombocompare_1_23s.out.tombocomparedone)
    tomboextract_23s_2 (tombocompare_2_23s.out.tombostat_ch, tombocompare_2_23s.out.tombocomparedone)
    tomboextract_23s_3 (tombocompare_3_23s.out.tombostat_ch, tombocompare_3_23s.out.tombocomparedone)


}