nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */


params.reference_16s = "/data/bhargava/reference_files/K12/k12_16s_88_extended.fa"
params.reference_23s = "/data/bhargava/reference_files/K12/k12_23s_78_extended.fa"

params.native16reads = "$baseDir/native/fastqs_16s_allmapped/bc1_subsample/*.fastq"
params.ivt16reads = "$baseDir/ivt/fastqs_16s_allmapped/bc1_subsample/*.fastq"

params.native23reads = "$baseDir/native/fastqs_23s_allmapped/bc1_subsample/*.fastq"
params.ivt23reads = "$baseDir/ivt/fastqs_23s_allmapped/bc1_subsample/*.fastq"

params.outdir_16s = "results/16s/bc1"
params.outdir_23s = "results/23s/bc1"

/*
 * Map fastqs to 16s rRNA reference using Minimap2 for native and IVT data
 */

process map_16s_native {

    publishDir "$params.outdir_16s/sams_16s_native", mode:'copy'
    
    tag "map fastq to get sam file"

    input:
    path reference_16s
    path nativefastqs

    output:
    path "*_16s.sam", emit: native16sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_16s $nativefastqs > ${nativefastqs.simpleName}_16s.sam --secondary=no
    """
}

process map_16s_ivt {

    publishDir "$params.outdir_16s/sams_16s_ivt", mode:'copy'

    tag "map fastq to get sam file"

    input:
    path reference_16s
    path ivtfastqs

    output:
    path "*_16s.sam", emit: ivt16sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_16s $ivtfastqs > ${ivtfastqs.simpleName}_16s.sam --secondary=no
    """
}

/*
 * Map fastqs to 23s rRNA reference using Minimap2 for native and IVT data
 */


process map_23s_native {

    publishDir "$params.outdir_23s/sams_23s_native", mode:'copy'

    tag "map fastq to get sam file"

    input:
    path reference_23s
    path nativefastqs

    output:
    path "*_23s.sam", emit: native23sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_23s $nativefastqs > ${nativefastqs.simpleName}_23s.sam --secondary=no
    """
}

process map_23s_ivt {

    publishDir "$params.outdir_23s/sams_23s_ivt", mode:'copy'
    
    tag "map fastq to get sam file"

    input:
    path reference_23s
    path ivtfastqs

    output:
    path "*_23s.sam", emit: ivt23sams

    script:

    """
    minimap2 -ax splice -uf -k14 $reference_23s $ivtfastqs > ${ivtfastqs.simpleName}_23s.sam --secondary=no
    """
}


/*
 * get mapping statistics for 16s sams for native and IVT data
 */

process flagstat_16s_native {

    publishDir "$params.outdir_16s/flagstat_16s_native", mode:'copy'
    
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

    publishDir "$params.outdir_16s/flagstat_16s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_23s/flagstat_23s_native", mode:'copy'
    
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

    publishDir "$params.outdir_23s/flagstat_23s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_16s/bams_16s_native", mode:'copy'
    
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

    publishDir "$params.outdir_16s/bams_16s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_23s/bams_23s_native", mode:'copy'
    
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

    publishDir "$params.outdir_23s/bams_23s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_16s/bams_16s_native", mode:'copy'
    
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

    publishDir "$params.outdir_16s/bams_16s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_23s/bams_23s_native", mode:'copy'
    
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

    publishDir "$params.outdir_23s/bams_23s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_16s/depths_16s_native", mode:'copy'
    
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

    publishDir "$params.outdir_16s/depths_16s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_23s/depths_23s_native", mode:'copy'
    
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

    publishDir "$params.outdir_23s/depths_23s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_16s/coverageplots16s_native", mode:'copy'
    
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

    publishDir "$params.outdir_16s/coverageplots16s_ivt", mode:'copy'
    
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

    publishDir "$params.outdir_23s/coverageplots23s_native", mode:'copy'
    
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

    publishDir "$params.outdir_23s/coverageplots23s_ivt", mode:'copy'
    
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
    .fromPath(params.native16reads)
    .set {native16fastqs_ch}

    Channel
    .fromPath(params.ivt16reads)
    .set {ivt16fastqs_ch}

    Channel
    .fromPath(params.native23reads)
    .set {native23fastqs_ch}

    Channel
    .fromPath(params.ivt23reads)
    .set {ivt23fastqs_ch}


    map_16s_native (reference_16s_ch, native16fastqs_ch)
    map_16s_ivt (reference_16s_ch, ivt16fastqs_ch)
    map_23s_native (reference_23s_ch, native23fastqs_ch)
    map_23s_ivt (reference_23s_ch, ivt23fastqs_ch)

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

}