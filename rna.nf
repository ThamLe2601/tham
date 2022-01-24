nextflow.enable.dsl=2
/* 
 * RNAseq thamle
 */
params.dir = "/home/thamle/rnaseq3"
params.outdir = "${params.dir}/output"
params.csvDir = '/home/thamle/rnaseq3/workflow/metadata.csv' 
params.fastqc = "${params.outdir}/fastqc"
params.trimed = "${params.outdir}/trim/"
params.genome = '/home/thamle/rnaseq3/ref/chr22_with_ERCC92.fa'
params.bed12 = '/home/thamle/rnaseq3/ref/chr22_with_ERCC92.bed12'
params.gtf = '/home/thamle/rnaseq3/ref/chr22_with_ERCC92.gtf'
params.index = "${params.outdir}/index"
params.star    = "${params.outdir}/star"
params.rseqc    = "${params.outdir}/rseqc"
params.qualimap    = "${params.outdir}/qualimap"
params.picard = "${params.outdir}/picard"
params.preseq = "${params.outdir}/preseq"
params.dupradar = "${params.outdir}/dupradar"
params.strandedness = 0
params.featurecount = "${params.outdir}/featurecount"
params.fc_extra_attributes = 'gene_name'
params.fc_group_features = 'gene_id'
params.fc_count_type = 'exon'
params.merge_featurecount = "${params.outdir}/merge_featurecount"
params.bam_suffix = "_Aligned.sortedByCoord.out.bam"
params.deseq2 = "${params.outdir}/deseq2"
params.design = '/home/thamle/rnaseq3/workflow/design.csv'
params.compare = '/home/thamle/rnaseq3/workflow/comparison.csv'
params.deseq2_fdr = 0.05
params.deseq2_lfc = 0.585
params.gprofiler = "${params.outdir}/gprofiler"
params.gprofiler_fdr = 0.05
params.multiqc = "${params.outdir}/multiqc"


log.info """\
 RNA__SEQ___NEXTFLOW

 """
  .stripIndent()
meta = Channel.fromPath(params.csvDir)
              .splitCsv(header:true)
              .map{ row-> tuple("$row.sampleid"), file("$row.read1"), file("$row.read2")}
              .set{sample_ch}
//genome_ch = Channel.fromPath(params.genome)
index_ch = Channel.fromPath(params.index)
                    .collect()
gtf_ch = Channel.fromPath(params.gtf)
                .collect()
bed12_ch = Channel.fromPath(params.bed12)
                  .collect()
design_csv = Channel.fromPath(params.design)
compare_ch = Channel.fromPath(params.compare)


process fastqc {
    tag "fatsqc $sampleid"
    publishDir "$params.fastqc", mode: 'copy'

    input:
    tuple val(sampleid), path(read1), path(read2)

    output:
    path '*_fastqc.{zip,html}', emit: report_qc
    path 'v_fastqc.txt', emit: version

    script:
    """
    fastqc --version &> v_fastqc.txt
    fastqc --quiet --threads $task.cpus $read1
    fastqc --quiet --threads $task.cpus $read2
    """
}

process trim_galore{
    tag "trim_galore $sampleid"
    publishDir "$params.trimed", mode: 'copy'

    input:
    tuple val(sampleid), path(read1), path(read2)

    output:
    tuple val(sampleid), path("*fq.gz"), emit: reads
    path '*trimming_report.txt', emit: report_trim
    path '*_fastqc.{zip,html}', emit: report_trim_qc

    
    
    //path "v_trim_galore.txt", emit: version//

    script:
    """
    trim_galore --paired --length 25 -q 30 --fastqc -j $task.cpus $read1 $read2 --gzip
    """
}

// process index {
//     tag "index"
//     publishDir"$params.index", mode: 'copy'

//     input:
//     path (genome_ch)
    
//     output:
//     file ("*")

//     script:
//     //mkdir index//
//     """
//     STAR --runThreadN 1 --genomeSAindexNbases 10 $task.cpus --runMode genomeGenerate --genomeFastaFiles $genome_ch --genomeDir ${params.index}
//     """
// }

process star {
    tag "star $sampleid"
    publishDir "$params.star", mode: 'copy'

    input:
    tuple val(sampleid), path(reads) 
    path ("${params.index}")


    output:
    tuple val(sampleid), path("${sampleid}_Aligned.sortedByCoord.out.bam"), emit: bam
    path("${sampleid}_Aligned.sortedByCoord.out.bam.{bai,csi}"), emit: bai
    path "*.out", emit: report_star
    path("*Log.final.out"), emit: log
    path "*SJ.out.tab"
    tuple val(meta), path("*Unmapped*"), emit: unmapped optional true
    //path "v_star.txt", emit: version//

    script:
    """ 
    STAR --readFilesIn $reads \
    --genomeDir ${params.index} \
    --outFileNamePrefix ${sampleid}_ \
    --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

    samtools index ${sampleid}_Aligned.sortedByCoord.out.bam
    """
}  

process rseqc {
    tag "rseqc $sampleid"
    publishDir "$params.rseqc", mode: 'copy'

    input:
    tuple val(sampleid), path(bam)
    path (bai)
    path (params.bed12)

    output:
    path "${sampleid}*", emit: report_rseqc
    //path "v_rseqc.txt", emit: version

    script:
    """
    read_distribution.py -i $bam -r ${params.bed12} > ${sampleid}.rseqc.read_distribution.txt
    read_duplication.py -i $bam -o ${sampleid}.rseqc.read_duplication
    inner_distance.py -i $bam -o ${sampleid}.rseqc -r ${params.bed12}
    infer_experiment.py -i $bam -r ${params.bed12} -s 2000000 > ${sampleid}.rseqc.infer_experiment.txt
    junction_saturation.py -i $bam -r ${params.bed12} -o ${sampleid}.rseqc
    bam_stat.py -i $bam > ${sampleid}.rseqc.bam_stat.txt
    junction_annotation.py -i $bam  -r ${params.bed12} -o ${sampleid} 2> ${sampleid}.rseqc.junction_annotation_log.txt
    """
}

process qualimap {
    tag "qualimap $sampleid"
    publishDir "$params.qualimap", mode: 'copy'

    input:
    tuple val(sampleid), path(bam)
    path(bai)
    path (params.gtf)

    output:
    path "${sampleid}*", emit: report_qualimap
    //path "v_qualimap.txt", emit: version

    script:
    """ 
    qualimap rnaseq -bam $bam -gtf ${params.gtf} -outfile ${sampleid}.pdf

    """
}

process picard {
    tag "picard $sampleid"
    publishDir "$params.picard", mode: 'copy'

    input:
    tuple val(sampleid), path (bam) 

    output:
    path ("${sampleid}*"), emit: report_picard

    script:
    """ 
    picard MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${sampleid}.markDups.bam \\
        METRICS_FILE=${sampleid}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

process preseq {
    tag "preseq $sampleid"
    publishDir "$params.preseq", mode: 'copy'

    input:
    tuple val(sampleid), path(bam)

    output:
    path ("${sampleid}*"), emit: report_preseq

    script:
    """ 
    preseq lc_extrap -v -B $bam -o ${sampleid}.ccurve.txt
    """
}

process dupradar {
    tag "dupradar $sampleid"
    publishDir "$params.dupradar", mode: 'copy'

    input:
    tuple val(sampleid), path(bam)
    path (params.gtf)

    output:
    path "*.{txt,pdf}", emit: report_dupradar
    // path "v_dupRadar.txt"

    script:

    """
    dupRadar.r $bam ${params.gtf} ${params.strandedness} paired ${task.cpus}
    Rscript -e "write(x=as.character(packageVersion('dupRadar')))"
    """

}

process featurecount {
    tag "featurecount $sampleid"
    publishDir "$params.featurecount", mode: 'copy'

    input:
    tuple val(sampleid), path (bam)
    path (params.gtf)

    output:
    path ("${sampleid}.featureCounts.txt"), emit: countout
    path "${sampleid}.featureCounts.txt.summary", emit: report_featurecount

    script:
    extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
    
    """ 
    featureCounts -a ${params.gtf} -g ${params.fc_group_features} -t ${params.fc_count_type} -s ${params.strandedness} -p -o ${sampleid}.featureCounts.txt $extraAttributes $bam 

    """
}

process merge_featurecount {
    tag "merge_featurecount"
    publishDir "$params.merge_featurecount", mode: 'copy'

    input:
    path (countout)

    output:
    path 'merged_gene_counts.txt', emit: merged_counts
    path 'gene_lengths.txt', emit: gene_lengths

    script:
    gene_ids = "<(tail -n +2 ${countout[0]} | cut -f1,7 )"
    counts_bam = countout.collect{filename ->
    // Remove first line and take 8th column (counts)
    "<(tail -n +2 ${filename} | sed 's:${params.bam_suffix}::' | cut -f8)"}.join(" ")

    """
    paste $gene_ids $counts_bam > merged_gene_counts.txt
    tail -n +2 ${countout[0]} | cut -f1,6 > gene_lengths.txt
    """
}

process deseq2 {
    tag "deseq2"
    publishDir "$params.deseq2", mode: 'copy'

    input:
    path(merged_counts)
    path (design_csv)
    path (compare_ch)
    
    output:
    path "*.{xlsx,jpg}", emit: download
    path "*_DESeq_results.tsv", emit: results_deseq2
    path "*{heatmap,plot,matrix}.tsv", emit: report_deseq2
    //path "v_DESeq2.txt", emit: version

    script:
    """ 
    DESeq2.r $merged_counts $design_csv $params.deseq2_fdr $params.deseq2_lfc $compare_ch
    Rscript -e "write(x=as.character(packageVersion('DESeq2')), file='v_DESeq2.txt')"
    """
}
process gprofiler {
    tag "deseq2"
    publishDir "$params.gprofiler", mode: 'copy'

    input:
    path (results_deseq2)

    output:
    path "*_gProfiler_results.tsv", emit: report_gprofiler
    path "*_gProfiler_results.xlsx", emit: download
    path "v_gProfiler.txt", emit: version

    script:
    """ 
    gProfiler.py $results_deseq2 -o hsapiens -q $params.deseq2_fdr -p $params.gprofiler_fdr
    pip freeze | grep gprofiler > v_gProfiler.txt
    """
}

process multiqc {
    tag "multiqc"
    publishDir "$params.multiqc", mode: 'copy'

    input:
    path report_qc
    path report_trim
    path report_trim_qc
    path report_star
    path report_rseqc
    path report_qualimap
    path report_picard
    path report_preseq
    path report_dupradar
    path report_featurecount
    path report_deseq2
    path report_gprofiler

    output:
    path 'multiqc_report.html'

    script:

    """ 
    multiqc .
    """
}

workflow {
   fastqc(sample_ch)
   trim_galore(sample_ch)
   //index(genome_ch)
   star(trim_galore.out.reads, index_ch) // index.out.collect(): error could not open genome file /home/thamle/rnaseq3/output/index/genomeParameters.txt
   rseqc(star.out.bam ,star.out.bai, bed12_ch)
   qualimap(star.out.bam ,star.out.bai, gtf_ch)
   picard(star.out.bam)
   preseq(star.out.bam)
   dupradar(star.out.bam, gtf_ch)
   featurecount(star.out.bam, gtf_ch)
   merge_featurecount(featurecount.out.countout.collect())
   deseq2(merge_featurecount.out.merged_counts, design_csv, compare_ch)
   gprofiler(deseq2.out.results_deseq2)
   multiqc(fastqc.out.report_qc.collect(), trim_galore.out.report_trim.collect(), trim_galore.out.report_trim_qc.collect(), star.out.report_star.collect(), rseqc.out.report_rseqc.collect(), qualimap.out.report_qualimap.collect(), picard.out.report_picard.collect(), preseq.out.report_preseq.collect(), dupradar.out.report_dupradar.collect(), featurecount.out.report_featurecount.collect(), deseq2.out.report_deseq2.collect(), gprofiler.out.report_gprofiler.collect())
}


