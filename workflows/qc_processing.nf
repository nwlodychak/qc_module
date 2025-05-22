#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { fastqc }   from '../modules/fastqc/main.nf'
include { posthoc_qc } from '../modules/post_hoc/main.nf'
/*
========================================================================================
    RUN WORKFLOWS
========================================================================================
*/

// Params
params.sample_sheet = '/Users/nwlodychak/venv/projects/qc_module/data/samplesheet.csv'
params.fastq_dir = '/Users/nwlodychak/venv/projects/qc_module/data'
params.qc_module_image_uri = "XXX"
params.qc_module_image_version = "0.1"
params.min_seqs = 1
// params.fastqc_mem = 2
// params.fastqc_cores = 2



workflow qc_processing {

    // Read samplesheet and create channel
    ch_samplesheet = Channel
        .fromPath("${params.sample_sheet}")
        .ifEmpty {
            throw new RuntimeException("Cannot find samplesheet! Check samplesheet path")
        }

    ch_fastq = Channel
        .fromPath("${params.fastq_dir}/*.fastq.gz")
        .map { fastq_file -> 
            def sample_id = fastq_file.simpleName
            tuple(sample_id, fastq_file)
    }

    fastqc( ch_fastq )
    println("COMPLETE!")
    // if (params.plotting) {
    posthoc_qc( fastqc.out.sample_id, fastqc.out.zip_file )
    // }
    // aggregate_qc()
}

workflow {
    qc_processing()
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}