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
 * get mapping statistics for 16s sams for native and IVT data
 */

process flagstat_16s_native {

    publishDir "$params.outdir/flagstat_16s_native", mode:'copy'
    
    tag "mapping stats with samtools flagstat"

    input:
    path native16sams

    output:
    path "*_flagstat.txt", emit: native16flagstat

    script:

    """
    samtools flagstat $native16sams | tee ${native16sams.simpleName}_flagstat.txt
    """
}

process flagstat_16s_ivt {

    publishDir "$params.outdir/flagstat_16s_ivt", mode:'copy'
    
    tag "mapping stats with samtools flagstat"

    input:
    path ivt16sams

    output:
    path "*_flagstat.txt", emit: ivt16flagstat

    script:

    """
    samtools flagstat $ivt16sams | tee ${ivt16sams.simpleName}_flagstat.txt
    """
}

/*
 * get mapping statistics for 23s sams for native and IVT data
 */

process flagstat_23s_native {

    publishDir "$params.outdir/flagstat_23s_native", mode:'copy'
    
    tag "mapping stats with samtools flagstat"

    input:
    path native23sams

    output:
    path "*_flagstat.txt", emit: native23flagstat

    script:

    """
    samtools flagstat $native23sams | tee ${native23sams.simpleName}_flagstat.txt
    """
}

process flagstat_23s_ivt {

    publishDir "$params.outdir/flagstat_23s_ivt", mode:'copy'
    
    tag "mapping stats with samtools flagstat"

    input:
    path ivt23sams

    output:
    path "*_flagstat.txt", emit: ivt23flagstat
    
    script:

    """
    samtools flagstat $ivt23sams | tee ${ivt23sams.simpleName}_flagstat.txt
    """
}

/*
 * Convert 16s sams into sorted bams for native and IVT data
 */

process bam_16s_native {

    publishDir "$params.outdir/bams_16s_native", mode:'copy'
    
    tag "convert sam to bam"

    input:
    path native16sams

    output:
    path "*_sorted.bam", emit: native16bams

    script:

    """
    samtools view -S -b -h $native16sams | samtools sort -o ${native16sams.simpleName}_sorted.bam
    """
}

process bam_16s_ivt {

    publishDir "$params.outdir/bams_16s_ivt", mode:'copy'
    
    tag "convert sam to bam"

    input:
    path ivt16sams

    output:
    path "*_sorted.bam", emit: ivt16bams

    script:

    """
    samtools view -S -b -h $ivt16sams | samtools sort -o ${ivt16sams.simpleName}_sorted.bam
    """
}

/*
 * Convert 23s sams into sorted bams for native and IVT data
 */

process bam_23s_native {

    publishDir "$params.outdir/bams_23s_native", mode:'copy'
    
    tag "convert sam to bam"

    input:
    path native23sams

    output:
    path "*_sorted.bam", emit: native23bams

    script:

    """
    samtools view -S -b -h $native23sams | samtools sort -o ${native23sams.simpleName}_sorted.bam
    """
}

process bam_23s_ivt {

    publishDir "$params.outdir/bams_23s_ivt", mode:'copy'
    
    tag "convert sam to bam"

    input:
    path ivt23sams

    output:
    path "*_sorted.bam", emit: ivt23bams

    script:

    """
    samtools view -S -b -h $ivt23sams | samtools sort -o ${ivt23sams.simpleName}_sorted.bam
    """
}

/*
 * Index 16s bams for native and IVT data
 */

process bamindex_16s_native {

    publishDir "$params.outdir/bams_16s_native", mode:'copy'
    
    tag "index the bam file"

    input:
    path native16bams

    output:
    path "*.bam.bai", emit: native16bamindex

    script:

    """
    samtools index $native16bams
    """
}

process bamindex_16s_ivt {

    publishDir "$params.outdir/bams_16s_ivt", mode:'copy'
    
    tag "index the bam file"

    input:
    path ivt16bams

    output:
    path "*.bam.bai", emit: ivt16bamindex

    script:

    """
    samtools index $ivt16bams
    """
}

/*
 * Index 23s bams for native and IVT data
 */

process bamindex_23s_native {

    publishDir "$params.outdir/bams_23s_native", mode:'copy'
    
    tag "index the bam file"

    input:
    path native23bams

    output:
    path "*.bam.bai", emit: native23bamindex

    script:

    """
    samtools index $native23bams
    """
}

process bamindex_23s_ivt {

    publishDir "$params.outdir/bams_23s_ivt", mode:'copy'
    
    tag "index the bam file"

    input:
    path ivt23bams

    output:
    path "*.bam.bai", emit: ivt23bamindex

    script:

    """
    samtools index $ivt23bams
    """
}

/*
 * Calculate depths from 16s bams for native and IVT data
 */

process depths_16s_native {

    publishDir "$params.outdir/depths_16s_native", mode:'copy'
    
    tag "calculate depth info from bam"

    input:
    path native16bams

    output:
    path "*.txt", emit: native16depths

    script:

    """
    samtools depth -a -m 0 $native16bams > ${native16bams.simpleName}.txt
    """
}

process depths_16s_ivt {

    publishDir "$params.outdir/depths_16s_ivt", mode:'copy'
    
    tag "calculate depth info from bam"

    input:
    path ivt16bams

    output:
    path "*.txt", emit: ivt16depths

    script:

    """
    samtools depth -a -m 0 $ivt16bams > ${ivt16bams.simpleName}.txt
    """
}

/*
 * Calculate depths from 23s bams for native and IVT data
 */


process depths_23s_native {

    publishDir "$params.outdir/depths_23s_native", mode:'copy'
    
    tag "calculate depth info from bam"

    input:
    path native23bams

    output:
    path "*.txt", emit: native23depths

    script:

    """
    samtools depth -a -m 0 $native23bams > ${native23bams.simpleName}.txt
    """
}

process depths_23s_ivt {

    publishDir "$params.outdir/depths_23s_ivt", mode:'copy'
    
    tag "calculate depth info from bam"

    input:
    path ivt23bams

    output:
    path "*.txt", emit: ivt23depths

    script:

    """
    samtools depth -a -m 0 $ivt23bams > ${ivt23bams.simpleName}.txt
    """
}

/*
 * Get coverage plots for 16s rRNA for native and IVT data
 */

process coverageplot_16s_native {

    publishDir "$params.outdir/coverageplots16s_native", mode:'copy'
    
    tag "coverage plot from depth file"

    input:
    path native16bams

    output:
    path "*.pdf", emit: native16coverageplot

    script:

    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    file = pd.read_csv("$native16bams", sep='\t', header=None)
    file.columns = ["chrom", "position", "depth"]
    sns.set_theme(style="darkgrid")
    x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
    plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
    plt.title("native - 16s rRNA")
    plt.xlabel("position on 16s rRNA")
    plt.ylabel("Coverage")
    plt.savefig("${native16bams.simpleName}" + ".pdf")
    """
}

process coverageplot_16s_ivt {

    publishDir "$params.outdir/coverageplots16s_ivt", mode:'copy'
    
    tag "coverage plot from depth file"

    input:
    path ivt16bams

    output:
    path "*.pdf", emit: ivt16coverageplot

    script:

    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    file = pd.read_csv("$ivt16bams", sep='\t', header=None)
    file.columns = ["chrom", "position", "depth"]
    sns.set_theme(style="darkgrid")
    x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
    plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
    plt.title("IVT - 16s rRNA")
    plt.xlabel("position on 16s rRNA")
    plt.ylabel("Coverage")
    plt.savefig("${ivt16bams.simpleName}" + ".pdf")
    """
}

/*
 * Get coverage plots for 23s rRNA for native and IVT data
 */

process coverageplot_23s_native {

    publishDir "$params.outdir/coverageplots23s_native", mode:'copy'
    
    tag "coverage plot from depth file"

    input:
    path native23bams

    output:
    path "*.pdf", emit: native23coverageplot

    script:

    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    file = pd.read_csv("$native23bams", sep='\t', header=None)
    file.columns = ["chrom", "position", "depth"]
    sns.set_theme(style="darkgrid")
    x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
    plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
    plt.title("native - 23s rRNA")
    plt.xlabel("position on 23s rRNA")
    plt.ylabel("Coverage")
    plt.savefig("${native23bams.simpleName}" + ".pdf")
    """
}

process coverageplot_23s_ivt {

    publishDir "$params.outdir/coverageplots23s_ivt", mode:'copy'
    
    tag "coverage plot from depth file"

    input:
    path ivt23bams

    output:
    path "*.pdf", emit: ivt23coverageplot

    script:

    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    file = pd.read_csv("$ivt23bams", sep='\t', header=None)
    file.columns = ["chrom", "position", "depth"]
    sns.set_theme(style="darkgrid")
    x = sns.relplot(x=file["position"], y=file["depth"], kind='line', height=6, aspect=4)
    plt.fill_between(file["position"], file["depth"], 0, facecolor="orange", color='blue', alpha=0.2)
    plt.title("IVT - 23s rRNA")
    plt.xlabel("position on 23s rRNA")
    plt.ylabel("Coverage")
    plt.savefig("${ivt23bams.simpleName}" + ".pdf")
    """
}

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

/*
 * Extract mapped reads into fastq files from sam files - 16s rRNA - native and IVT data using samtools
 */

process extract_16s_native {

    publishDir "$params.outdir/mappedfastqs_16s_native", mode:'copy'
    
    tag "extract mapped reads from sam files into fastq"

    input:
    path native16sams

    output:
    path "*.fastq", emit: native16mappedfastqs

    script:

    """
    samtools fastq -F4 $native16sams > ${native16sams.simpleName}.fastq
    """
}

process extract_16s_ivt {

    publishDir "$params.outdir/mappedfastqs_16s_ivt", mode:'copy'
    
    tag "extract mapped reads from sam files into fastq"

    input:
    path ivt16sams

    output:
    path "*.fastq", emit: ivt16mappedfastqs

    script:

    """
    samtools fastq -F4 $ivt16sams > ${ivt16sams.simpleName}.fastq
    """
}

/*
 * Extract mapped reads into fastq files from sam files - 23s rRNA - native and IVT data using samtools
 */

process extract_23s_native {

    publishDir "$params.outdir/mappedfastqs_23s_native", mode:'copy'
    
    tag "extract mapped reads from sam files into fastq"

    input:
    path native23sams

    output:
    path "*.fastq", emit: native23mappedfastqs

    script:

    """
    samtools fastq -F4 $native23sams > ${native23sams.simpleName}.fastq
    """
}

process extract_23s_ivt {

    publishDir "$params.outdir/mappedfastqs_23s_ivt", mode:'copy'
    
    tag "extract mapped reads from sam files into fastq"

    input:
    path ivt23sams

    output:
    path "*.fastq", emit: ivt23mappedfastqs

    script:

    """
    samtools fastq -F4 $ivt23sams > ${ivt23sams.simpleName}.fastq
    """
}

/*
 * Mapped native 16s fastqs to sorted bams
 */

process mappedbams_16s_native {

    publishDir "$params.outdir/mappedbams_16s_native", mode:'copy'
    
    tag "mapped 16s native fastqs to sorted bam files"

    input:
    path reference_16s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: native16mappedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_16s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}

/*
 * Mapped native 23s fastqs to sorted bams
 */

process mappedbams_23s_native {

    publishDir "$params.outdir/mappedbams_23s_native", mode:'copy'
    
    tag "mapped 23s native fastqs to sorted bam files"

    input:
    path reference_23s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: native23mappedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_23s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}

/*
 * Mapped IVT 16s fastqs to sorted bams
 */

process mappedbams_16s_ivt {

    publishDir "$params.outdir/mappedbams_16s_ivt", mode:'copy'
    
    tag "mapped 16s IVT fastqs to sorted bam files"

    input:
    path reference_16s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: ivt16mappedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_16s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}

/*
 * Mapped IVT 23s fastqs to sorted bams
 */

process mappedbams_23s_ivt {

    publishDir "$params.outdir/mappedbams_23s_ivt", mode:'copy'
    
    tag "mapped 23s IVT fastqs to sorted bam files"

    input:
    path reference_23s
    path mappedfastqs

    output:
    path "*_sorted.bam", emit: ivt23mappedbams

    script:

    """
    minimap2 -ax splice -uf -k14 --secondary=no $reference_23s $mappedfastqs | samtools view -S -b -h | samtools sort -o ${mappedfastqs.simpleName}_sorted.bam
    """
}


/*
 * Index mapped 16s bams for native and IVT data
 */

process mappedbamindex_16s_native {

    publishDir "$params.outdir/mappedbams_16s_native", mode:'copy'
    
    tag "index the mapped bam file"

    input:
    path native16bams

    output:
    path "*.bam.bai", emit: native16mappedbamindex

    script:

    """
    samtools index $native16bams
    """
}

process mappedbamindex_16s_ivt {

    publishDir "$params.outdir/mappedbams_16s_ivt", mode:'copy'
    
    tag "index the mapped bam file"

    input:
    path ivt16bams

    output:
    path "*.bam.bai", emit: ivt16mappedbamindex

    script:

    """
    samtools index $ivt16bams
    """
}

/*
 * Index mapped 23s bams for native and IVT data
 */

process mappedbamindex_23s_native {

    publishDir "$params.outdir/mappedbams_23s_native", mode:'copy'
    
    tag "index the mapped bam file"

    input:
    path native23bams

    output:
    path "*.bam.bai", emit: native23mappedbamindex

    script:

    """
    samtools index $native23bams
    """
}

process mappedbamindex_23s_ivt {

    publishDir "$params.outdir/mappedbams_23s_ivt", mode:'copy'
    
    tag "index the mapped bam file"

    input:
    path ivt23bams

    output:
    path "*.bam.bai", emit: ivt23mappedbamindex

    script:

    """
    samtools index $ivt23bams
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