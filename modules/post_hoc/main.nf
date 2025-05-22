#!/usr/bin/env nextflow

// Process
process posthoc_qc {
    tag "Post-Hoc Processing for ${sample_id}..."

    input:
    val(sample_id)
    path(zip_file)
    
    output:
    val(sample_id), emit: sample_id
    path('plots/*.png'), emit: plots
    path('tables/*.csv'), emit: tables

    errorStrategy 'retry'
    maxRetries 1
    publishDir "${workflow.launchDir}/results/by_sample_results/${sample_id}/posthoc/", mode: 'copy', failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.posthoc_cores ?: 2
    memory params.posthoc_mem ?: '2.GB'
    
    script:
    """
    Rscript bin/qc_module ${sample_id} ${zip_file} ${params.samplesheet} ${params.min_seqs}
    """
}

workflow {
    posthoc_qc(fqc_ch)
}

