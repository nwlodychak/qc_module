#!/usr/bin/env nextflow

// Process
process fastqc {
    tag "FastQC for ${sample_id}"
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${workflow.launchDir}/results/by_sample_results/${sample_id}/fastqc/", mode: 'copy', failOnError: true
    container "${params.qc_module_image_uri}:${params.qc_module_image_version}"
    cpus params.fastqc_cores ?: 2
    memory params.fastqc_mem ?: '2.GB'

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    val(sample_id), emit: sample_id
    tuple val(sample_id), path("*.zip"), emit: zip_file
    tuple val(sample_id), path("*.html"), emit: html_file
    
    script:
    """
    fastqc ${fastq} -o .
    """
}