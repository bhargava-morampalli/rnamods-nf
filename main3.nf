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
 * Import modules
 */



include { nanodocpreprocess as nanodocpreprocess_n_1_16s; nanodocpreprocess as nanodocpreprocess_n_2_16s; nanodocpreprocess as nanodocpreprocess_n_3_16s; nanodocpreprocess as nanodocpreprocess_n_1_23s; nanodocpreprocess as nanodocpreprocess_n_2_23s; nanodocpreprocess as nanodocpreprocess_n_3_23s } from '/home/bhargavam/Documents/nextflowmodules/nanodocpreprocess'

include { nanodocpreprocess as nanodocpreprocess_i_1_16s; nanodocpreprocess as nanodocpreprocess_i_2_16s; nanodocpreprocess as nanodocpreprocess_i_3_16s; nanodocpreprocess as nanodocpreprocess_i_1_23s; nanodocpreprocess as nanodocpreprocess_i_2_23s; nanodocpreprocess as nanodocpreprocess_i_3_23s } from '/home/bhargavam/Documents/nextflowmodules/nanodocpreprocess'

include { nanodoccompare_16s as nanodoccompare_1_16s; nanodoccompare_16s as nanodoccompare_2_16s; nanodoccompare_16s as nanodoccompare_3_16s } from '/home/bhargavam/Documents/nextflowmodules/nanodoccompare_16s'

include { nanodoccompare_23s as nanodoccompare_1_23s; nanodoccompare_23s as nanodoccompare_2_23s; nanodoccompare_23s as nanodoccompare_3_23s } from '/home/bhargavam/Documents/nextflowmodules/nanodoccompare_23s'

include { xporeprep as xporeprep1_n_16s; xporeprep as xporeprep2_n_16s; xporeprep as xporeprep3_n_16s } from '/home/bhargavam/Documents/nextflowmodules/xporeprep'

include { xporeprep as xporeprep1_n_23s; xporeprep as xporeprep2_n_23s; xporeprep as xporeprep3_n_23s } from '/home/bhargavam/Documents/nextflowmodules/xporeprep'

include { xporeprep as xporeprep1_i_16s; xporeprep as xporeprep2_i_16s; xporeprep as xporeprep3_i_16s } from '/home/bhargavam/Documents/nextflowmodules/xporeprep'

include { xporeprep as xporeprep1_i_23s; xporeprep as xporeprep2_i_23s; xporeprep as xporeprep3_i_23s } from '/home/bhargavam/Documents/nextflowmodules/xporeprep'

include { yanocompprep as yanocompprep1_n_16s; yanocompprep as yanocompprep2_n_16s; yanocompprep as yanocompprep3_n_16s } from '/home/bhargavam/Documents/nextflowmodules/yanocompprep'

include { yanocompprep as yanocompprep1_n_23s; yanocompprep as yanocompprep2_n_23s; yanocompprep as yanocompprep3_n_23s } from '/home/bhargavam/Documents/nextflowmodules/yanocompprep'

include { yanocompprep as yanocompprep1_i_16s; yanocompprep as yanocompprep2_i_16s; yanocompprep as yanocompprep3_i_16s } from '/home/bhargavam/Documents/nextflowmodules/yanocompprep'

include { yanocompprep as yanocompprep1_i_23s; yanocompprep as yanocompprep2_i_23s; yanocompprep as yanocompprep3_i_23s } from '/home/bhargavam/Documents/nextflowmodules/yanocompprep'

include { yanocompanalysis as yanocompanalysis_1_16s; yanocompanalysis as yanocompanalysis_2_16s; yanocompanalysis as yanocompanalysis_3_16s } from '/home/bhargavam/Documents/nextflowmodules/yanocompanalysis'

include { yanocompanalysis as yanocompanalysis_1_23s; yanocompanalysis as yanocompanalysis_2_23s; yanocompanalysis as yanocompanalysis_3_23s } from '/home/bhargavam/Documents/nextflowmodules/yanocompanalysis'


include { xporeanalysis as xporeanalysis_1_16s; xporeanalysis as xporeanalysis_2_16s; xporeanalysis as xporeanalysis_3_16s } from '/home/bhargavam/Documents/nextflowmodules/xporeanalysis'

include { xporeanalysis as xporeanalysis_1_23s; xporeanalysis as xporeanalysis_2_23s; xporeanalysis as xporeanalysis_3_23s } from '/home/bhargavam/Documents/nextflowmodules/xporeanalysis'


/*
 * WORKFLOW
 */


workflow {


    Channel
    .value(file(params.reference_16s))
    .set {reference_16s_ch}

    Channel
    .value(file(params.reference_23s))
    .set {reference_23s_ch}


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
    .value(file("$baseDir/results/fast5s_native_1_16s_single_ndocpreprocess/"))
    .set {native_1_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_native_2_16s_single_ndocpreprocess/"))
    .set {native_2_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_native_3_16s_single_ndocpreprocess/"))
    .set {native_3_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_native_1_23s_single_ndocpreprocess/"))
    .set {native_1_23s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_native_2_23s_single_ndocpreprocess/"))
    .set {native_2_23s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_native_3_23s_single_ndocpreprocess/"))
    .set {native_3_23s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_1_16s_single_ndocpreprocess/"))
    .set {ivt_1_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_2_16s_single_ndocpreprocess/"))
    .set {ivt_2_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_3_16s_single_ndocpreprocess/"))
    .set {ivt_3_16s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_1_23s_single_ndocpreprocess/"))
    .set {ivt_1_23s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_2_23s_single_ndocpreprocess/"))
    .set {ivt_2_23s_nanodoc}

    Channel
    .value(file("$baseDir/results/fast5s_ivt_3_23s_single_ndocpreprocess/"))
    .set {ivt_3_23s_nanodoc}








    Channel
    .value(file("$baseDir/results/native_1_16s_f5c/*1*eventalign.txt"))
    .set {native_1_16s_eventalign}

    Channel
    .value(file("$baseDir/results/native_2_16s_f5c/*2*eventalign.txt"))
    .set {native_2_16s_eventalign}

    Channel
    .value(file("$baseDir/results/native_3_16s_f5c/*3*eventalign.txt"))
    .set {native_3_16s_eventalign}

    Channel
    .value(file("$baseDir/results/native_1_23s_f5c/*1*eventalign.txt"))
    .set {native_1_23s_eventalign}

    Channel
    .value(file("$baseDir/results/native_2_23s_f5c/*2*eventalign.txt"))
    .set {native_2_23s_eventalign}

    Channel
    .value(file("$baseDir/results/native_3_23s_f5c/*3*eventalign.txt"))
    .set {native_3_23s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_1_16s_f5c/*1*eventalign.txt"))
    .set {ivt_1_16s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_2_16s_f5c/*2*eventalign.txt"))
    .set {ivt_2_16s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_3_16s_f5c/*3*eventalign.txt"))
    .set {ivt_3_16s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_1_23s_f5c/*1*eventalign.txt"))
    .set {ivt_1_23s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_2_23s_f5c/*2*eventalign.txt"))
    .set {ivt_2_23s_eventalign}

    Channel
    .value(file("$baseDir/results/ivt_3_23s_f5c/*3*eventalign.txt"))
    .set {ivt_3_23s_eventalign}



    Channel
    .value(file("$baseDir/results/native_1_16s_f5c/*1*summary.txt"))
    .set {native_1_16s_summary}

    Channel
    .value(file("$baseDir/results/native_2_16s_f5c/*2*summary.txt"))
    .set {native_2_16s_summary}

    Channel
    .value(file("$baseDir/results/native_3_16s_f5c/*3*summary.txt"))
    .set {native_3_16s_summary}

    Channel
    .value(file("$baseDir/results/native_1_23s_f5c/*1*summary.txt"))
    .set {native_1_23s_summary}

    Channel
    .value(file("$baseDir/results/native_2_23s_f5c/*2*summary.txt"))
    .set {native_2_23s_summary}

    Channel
    .value(file("$baseDir/results/native_3_23s_f5c/*3*summary.txt"))
    .set {native_3_23s_summary}

    Channel
    .value(file("$baseDir/results/ivt_1_16s_f5c/*1*summary.txt"))
    .set {ivt_1_16s_summary}

    Channel
    .value(file("$baseDir/results/ivt_2_16s_f5c/*2*summary.txt"))
    .set {ivt_2_16s_summary}

    Channel
    .value(file("$baseDir/results/ivt_3_16s_f5c/*3*summary.txt"))
    .set {ivt_3_16s_summary}

    Channel
    .value(file("$baseDir/results/ivt_1_23s_f5c/*1*summary.txt"))
    .set {ivt_1_23s_summary}

    Channel
    .value(file("$baseDir/results/ivt_2_23s_f5c/*2*summary.txt"))
    .set {ivt_2_23s_summary}

    Channel
    .value(file("$baseDir/results/ivt_3_23s_f5c/*3*summary.txt"))
    .set {ivt_3_23s_summary}




    Channel
    .value(file("$baseDir/data/16s_yano_1.yml"))
    .set {yml_16s_1}

    Channel
    .value(file("$baseDir/data/16s_yano_2.yml"))
    .set {yml_16s_2}

    Channel
    .value(file("$baseDir/data/16s_yano_3.yml"))
    .set {yml_16s_3}

    Channel
    .value(file("$baseDir/data/23s_yano_1.yml"))
    .set {yml_23s_1}

    Channel
    .value(file("$baseDir/data/23s_yano_2.yml"))
    .set {yml_23s_2}

    Channel
    .value(file("$baseDir/data/23s_yano_3.yml"))
    .set {yml_23s_3}




    nanodocpreprocess_n_1_16s (native_1_16s_single_fast5s_ch, reference_16s_ch)
    nanodocpreprocess_n_2_16s (native_2_16s_single_fast5s_ch, reference_16s_ch)
    nanodocpreprocess_n_3_16s (native_3_16s_single_fast5s_ch, reference_16s_ch)

    nanodocpreprocess_n_1_23s (native_1_23s_single_fast5s_ch, reference_23s_ch)
    nanodocpreprocess_n_2_23s (native_2_23s_single_fast5s_ch, reference_23s_ch)
    nanodocpreprocess_n_3_23s (native_3_23s_single_fast5s_ch, reference_23s_ch)

    nanodocpreprocess_i_1_16s (ivt_1_16s_single_fast5s_ch, reference_16s_ch)
    nanodocpreprocess_i_2_16s (ivt_2_16s_single_fast5s_ch, reference_16s_ch)
    nanodocpreprocess_i_3_16s (ivt_3_16s_single_fast5s_ch, reference_16s_ch)

    nanodocpreprocess_i_1_23s (ivt_1_23s_single_fast5s_ch, reference_23s_ch)
    nanodocpreprocess_i_2_23s (ivt_2_23s_single_fast5s_ch, reference_23s_ch)
    nanodocpreprocess_i_3_23s (ivt_3_23s_single_fast5s_ch, reference_23s_ch)

    yanocompprep1_n_16s (native_1_16s_eventalign, native_1_16s_summary)
    yanocompprep2_n_16s (native_2_16s_eventalign, native_2_16s_summary)
    yanocompprep3_n_16s (native_3_16s_eventalign, native_3_16s_summary)

    yanocompprep1_i_16s (ivt_1_16s_eventalign, ivt_1_16s_summary)
    yanocompprep2_i_16s (ivt_2_16s_eventalign, ivt_2_16s_summary)
    yanocompprep3_i_16s (ivt_3_16s_eventalign, ivt_3_16s_summary)

    yanocompprep1_n_23s (native_1_23s_eventalign, native_1_23s_summary)
    yanocompprep2_n_23s (native_2_23s_eventalign, native_2_23s_summary)
    yanocompprep3_n_23s (native_3_23s_eventalign, native_3_23s_summary)

    yanocompprep1_i_23s (ivt_1_23s_eventalign, ivt_1_23s_summary)
    yanocompprep2_i_23s (ivt_2_23s_eventalign, ivt_2_23s_summary)
    yanocompprep3_i_23s (ivt_3_23s_eventalign, ivt_3_23s_summary)

    yanocompanalysis_1_16s (yanocompprep1_n_16s.out.yanohdf5, yanocompprep1_i_16s.out.yanohdf5, yanocompprep1_n_16s.out.yanoprepdone, yanocompprep1_i_16s.out.yanoprepdone)
    yanocompanalysis_2_16s (yanocompprep2_n_16s.out.yanohdf5, yanocompprep2_i_16s.out.yanohdf5, yanocompprep2_n_16s.out.yanoprepdone, yanocompprep2_i_16s.out.yanoprepdone)
    yanocompanalysis_3_16s (yanocompprep3_n_16s.out.yanohdf5, yanocompprep3_i_16s.out.yanohdf5, yanocompprep3_n_16s.out.yanoprepdone, yanocompprep3_i_16s.out.yanoprepdone)

    yanocompanalysis_1_23s (yanocompprep1_n_23s.out.yanohdf5, yanocompprep1_i_23s.out.yanohdf5, yanocompprep1_n_23s.out.yanoprepdone, yanocompprep1_i_23s.out.yanoprepdone)
    yanocompanalysis_2_23s (yanocompprep2_n_23s.out.yanohdf5, yanocompprep2_i_23s.out.yanohdf5, yanocompprep2_n_23s.out.yanoprepdone, yanocompprep2_i_23s.out.yanoprepdone)
    yanocompanalysis_3_23s (yanocompprep3_n_23s.out.yanohdf5, yanocompprep3_i_23s.out.yanohdf5, yanocompprep3_n_23s.out.yanoprepdone, yanocompprep3_i_23s.out.yanoprepdone)

    xporeprep1_n_16s (native_1_16s_eventalign)
    xporeprep2_n_16s (native_2_16s_eventalign)
    xporeprep3_n_16s (native_3_16s_eventalign)

    xporeprep1_i_16s (ivt_1_16s_eventalign)
    xporeprep2_i_16s (ivt_2_16s_eventalign)
    xporeprep3_i_16s (ivt_3_16s_eventalign)

    xporeprep1_n_23s (native_1_23s_eventalign)
    xporeprep2_n_23s (native_2_23s_eventalign)
    xporeprep3_n_23s (native_3_23s_eventalign)

    xporeprep1_i_23s (ivt_1_23s_eventalign)
    xporeprep2_i_23s (ivt_2_23s_eventalign)
    xporeprep3_i_23s (ivt_3_23s_eventalign)

    xporeanalysis_1_16s (yml_16s_1, xporeprep1_n_16s.out.xporeprepdone, xporeprep1_i_16s.out.xporeprepdone)
    xporeanalysis_2_16s (yml_16s_2, xporeprep2_n_16s.out.xporeprepdone, xporeprep2_i_16s.out.xporeprepdone)
    xporeanalysis_3_16s (yml_16s_3, xporeprep3_n_16s.out.xporeprepdone, xporeprep3_i_16s.out.xporeprepdone)

    xporeanalysis_1_23s (yml_23s_1, xporeprep1_n_23s.out.xporeprepdone, xporeprep1_i_23s.out.xporeprepdone)
    xporeanalysis_2_23s (yml_23s_2, xporeprep2_n_23s.out.xporeprepdone, xporeprep2_i_23s.out.xporeprepdone)
    xporeanalysis_3_23s (yml_23s_3, xporeprep3_n_23s.out.xporeprepdone, xporeprep3_i_23s.out.xporeprepdone)

}