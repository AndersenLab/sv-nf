#!/usr/bin/env nextflow

// Use DSL2
nextflow.preview.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

/*
~ ~ ~ > * PARAMETERS SETUP
*/

// Variables
date = new Date().format('yyyyMMdd')

// Setup pipeline parameter


/*
~ ~ ~ > * LOG AND HELP MESSAGE SETUP
*/

if (!params.help) {
log.info '''
SV - N F   P I P E L I N E
===============================================
'''
    log.info ""
    log.info "Project           = ${params.project}"
    log.info "CP pipeline       = ${params.pipeline}"
    log.info "Groups            = ${params.groups}"
    log.info "Output            = ${params.out}"
    log.info ""
    } else {
log.info '''
SV - N F   P I P E L I N E
===============================================
'''
    log.info "Usage:"
    log.info "The typical command for running the pipeline is as follows:"
    log.info "nextflow run main.nf --pipeline <CellProfiler pipeline to use> --project <path to your project directory>"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--project                      The path to your project directory"
    log.info "--pipeline                     The CP pipeline to use: toxin, dauer"
    log.info ""
    log.info "Optional arguments:"
    log.info "--groups                       comma separated metadata groupings for CellProfiler, default is plate,well"
    log.info "--outdir                       Output directory to place files, default is project/Analysis-{current date}"
    log.info "--help                         This usage statement."
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/

workflow {
    \\ add stuff here
}

/*
~ ~ ~ > * CONFIGURE FILES FOR CELLPROFILER
*/

process config_CP_input_dauer {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    input:
        tuple file(raw_pipe), val(meta_dir), val(meta), val(model_dir), val(model1), val(model2),
        file(config_script), val(project), val(mask), val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}

process config_CP_input_toxin {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    input:
        tuple file(raw_pipe), val(meta_dir), val(meta), val(model_dir), val(model1), val(model2), val(model3), val(model4),
        file(config_script), val(project), val(mask), val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' | \\
        awk '{gsub(/MODEL3_XML_FILE/,"${model3}"); print}' | \\
        awk '{gsub(/MODEL4_XML_FILE/,"${model4}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}

/*
~ ~ ~ > * RUN CELLPROFILER
*/

process runCP {

    label "cellpro"

    input:
        tuple val(group), file(pipeline), file(output)

    output:
        stdout emit: cp_output //tuple file("*.csv"), file("*.png"), emit: cp_output

    """
        # Run cellprofiler headless
        cellprofiler -c -r -p ${pipeline} \
        -g ${group} \
        -o ${output}

    """
}

/*
~ ~ ~ > * PROCESS CELLPROFILER OUTPUTS
*/

process proc_CP_output_dauer {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
        tuple val(cp_output), val(out_dir), val(model_name1), val(model_name2), file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data

        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images

        # find .csv files, concatenate them, and write new file
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name1}.csv > ${out_dir}/processed_data/${model_name1}.csv
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name2}.csv > ${out_dir}/processed_data/${model_name2}.csv
        
        # move all the output images to process_images directory
        mv ${out_dir}/CP_output/**/*.png ${out_dir}/processed_images

        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}

process proc_CP_output_toxin {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
        tuple val(cp_output), val(out_dir), val(model_name1), val(model_name2),
        val(model_name3), val(model_name4), file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data

        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images

        # find .csv files, concatenate them, and write new file
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name1}.csv > ${out_dir}/processed_data/${model_name1}.csv
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name2}.csv > ${out_dir}/processed_data/${model_name2}.csv
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name3}.csv > ${out_dir}/processed_data/${model_name3}.csv
        awk 'FNR==1 && NR!=1 { while (/^ImageNumber/) getline; } 1 {print}' ${out_dir}/CP_output/*/${model_name4}.csv > ${out_dir}/processed_data/${model_name4}.csv
        
        # move all the output images to process_images directory
        mv ${out_dir}/CP_output/**/*.png ${out_dir}/processed_images

        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}

/*
~ ~ ~ > * GENERATE REPORT
*/
workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    { Parameters }
    ---------------------------
    Project                                 = ${params.project}
    Pipeline Used                           = ${params.pipeline}
    Result Directory                        = ${params.out}
    """

    println summary

}