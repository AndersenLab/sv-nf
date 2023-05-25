#!/usr/bin/env nextflow

// Use DSL2
nextflow.preview.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

// Variables
date = new Date().format('yyyyMMdd')

// Parameters
params.release = null
params.help = null
params.debug = null
params.bin_dir = "${workflow.projectDir}/bin"

// Parameter logic
if(params.debug) {
    println """
        *** Using debug mode ***
    """
    params.species = "elegans"
    params.bam_dir = "${workflow.projectDir}/debug" // path to directory holding .bam files
    params.out = "SV_results_DEBUG_${params.species}_${date}"
    params.ref = "${workflow.projectDir}/debug/debug_ref.fa.gz" 
    params.sp_sheet = "${workflow.projectDir}/debug/sample_sheet.txt"
} else {
    params.species = null
    params.bam_dir = "/projects/b1059/data/c_${params.species}/WI/alignments" // path to directory holding .bam files
    params.out = "SV_results_${params.species}_${date}"
    params.ref = null
    params.sp_sheet = null 
}

// LOG AND HELP MESSAGE SETUP
if (!params.help) {
log.info '''
S V - N F    P I P E L I N E
===============================================
'''
    log.info ""
    log.info "Species                  = ${params.species}"
    log.info "Reference File           = ${params.ref}"
    log.info "CaeNDR Release           = ${params.release}"
    log.info "Sample Sheet             = ${params.sp_sheet}"
    log.info ".bam Directory           = ${params.bam_dir}"
    log.info "Output Directory         = ${params.out}"
    log.info "Debug                    = ${params.debug}"
    log.info ""
    } else {
log.info '''
S V - N F    P I P E L I N E
===============================================
'''
    log.info "Usage:"
    log.info "The typical command for running the pipeline is as follows:"
    log.info "nextflow run main.nf --species elegans --ref /projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa --release <latest CaeNDR release>"
    log.info "To debug use:"
    log.info "nextflow run main.nf --debug"
    log.info ""
    log.info "Required Arguments:"
    log.info "--species     String           One of three, elegans, briggsae, tropicalis"
    log.info "--ref         String           Full path to the .fa uncompressed reference file, default set to null"
    log.info "--release     String           The 8-digit date code for CaeNDR release, e.g 20220216. Required if --sp_sheet not specified"
    log.info "OR"
    log.info "--sp_sheet    String           A path to the sample_sheet.txt file for calling INDELs instead of release. Required if --release not specified"
    log.info "--PIF         BOOL            To perfrom pairwise-indel finder (PIF) analysis, default is false"
    log.info ""
    log.info "Optional Arguments:"
    log.info "--bam_dir     String           The path to the .bam directory, default set for QUEST: /projects/b1059/data/c_<species>/WI/alignments"
    log.info "--out         String           The output directory, default is SV_indel_results_<species>_<date>"
    log.info "--debug"
    log.info ""
    log.info "Flags:"
    log.info "--help                          Display this message"
    log.info ""
    log.info "Notes:"
    log.info "The briggsae reference path on QUEST is /projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa.gz"
    log.info "The tropicalis reference path on QUEST is /projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa.gz"
    log.info "--------------------------------------------------------"
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/
if ( params.PIF==null ) println "Using DELLY all mode" 

workflow {

    // Choose how we get isotype reference strains
    if(params.release){

        delly_in = Channel.fromPath("${params.bin_dir}/config_delly.R")
            .combine(Channel.from("${params.species}")) // get strain names from WI sheets
            .combine(Channel.from("${params.release}")) // get strain names from WI sheets
            .combine(Channel.from("${params.bam_dir}"))
            .combine(Channel.from("${params.ref}"))
            //.view()
    } else {
        delly_in = Channel.fromPath("${params.bin_dir}/config_delly.R")
            .combine(Channel.from("${params.species}")) // get strain names from WI sheets
            .combine(Channel.from("${params.sp_sheet}")) // take strain names from sample sheet
            .combine(Channel.from("${params.bam_dir}"))
            .combine(Channel.from("${params.ref}"))
            //.view()
    }


    // make the delly run parameters
    config_delly(delly_in)

    // setup the channel to run delly
    delly_pif_ch = config_delly.out.delly_par_file
        .splitCsv(header:true, sep: "\t")
            .map { row -> [row.strain, file("${row.bam}"), file("${row.index}"), file("${row.ref}")] }
            //.view()

    // run delly for pairwise-indel finder (pif)
    
    if ( params.PIF==null | params.PIF==false ){
        
        // run delly calling all SV types
        delly_all(delly_pif_ch)
        // collect output files 
        merge_delly_all_ch = delly_all.out.collect()
        
        // merge delly output files
        merge_delly_all(merge_delly_all_ch)

        // Collect outputs for genotype
        geno_sites_ch = delly_pif_ch.combine(merge_delly_all.out)


    } else {
        // run delly for Pairwise-indel finder (PIF)
        delly_pif(delly_pif_ch)
        // collect output files
        merge_delly_pif_ch = delly_pif.out.collect()

        // merge delly output files
        merge_delly_pif(merge_delly_pif_ch)

        // Collect outputs for genotype
        geno_sites_ch = delly_pif_ch.combine(merge_delly_pif.out)

    }

    // run genotype_sites
    genotype_sites(geno_sites_ch)

    // send output to proccess delly channel
    if ( params.PIF==null | params.PIF==false ){
        proc_genos_ch = genotype_sites.out.delly_genos.collect() | proc_genos_all
    
    } else {
    
        proc_genos_ch = genotype_sites.out.delly_genos.collect() | proc_genos
            // setup ceandr_pif script and inputs as a channel
        ceandr_pif_ch = Channel.fromPath("${params.bin_dir}/bed_to_VCF.R")
        .combine(proc_genos.out.raw_ceandr_bed)
        .combine(Channel.fromPath("${params.ref}"))
        .combine(Channel.from("${params.species}"))

        // run it
        output_caendr_pif(ceandr_pif_ch)
    }
        //.view()


}

process config_delly {
    
    label "R"

    input:
        tuple file(config_script), val(species), val(samples), val(bams_path), val(ref_path)

    output:
        path "sample_file.txt", emit: sample_file
        path "delly_config_file.tsv", emit: delly_par_file


    """
        # Use config script to setup delly run parameters
        Rscript --vanilla ${config_script} ${species} ${samples} ${bams_path} ${ref_path}

    """
}

process delly_pif {
    
    label "dell_big"

    input:
        tuple val(strain), file(bam), file(index), file(ref)

    output:
        file "*.bcf"


    """
        delly call -t DEL ${bam} -g ${ref} -o ${strain}_del.bcf
        delly call -t INS ${bam} -g ${ref} -o ${strain}_ins.bcf
        
    """
}

process delly_all {
    
    label "dell_big"

    input:
        tuple val(strain), file(bam), file(index), file(ref)

    output:
        file "*.bcf"


    """
        delly call ${bam} -g ${ref} -o ${strain}_.bcf
    """
}

process merge_delly_pif {
    
    label "dell"

    input:
        path("*")

    output:
        file "*.bcf" 

    """
        vcf_list=`echo *.bcf`
        delly merge \${vcf_list} --minsize 50 --maxsize 500 -b 1000 -r 0.8000012 -o WI.merge.sites.bcf
    """

}

process merge_delly_all {
    
    label "dell"

    input:
        path("*")

    output:
        file "*.bcf" 

    """
        vcf_list=`echo *.bcf`
        delly merge \${vcf_list} -b 1000 -r 0.8000012 -o WI.merge.sites.bcf
    """

}

process genotype_sites {

    label "dell_big"
        
    input:
        tuple val(isotype), file(bam), file(index), file(ref), file(merged_bcf)

    output:
        tuple file("${isotype}_geno.bcf"), file("${isotype}_geno.bcf.csi"), emit: delly_genos

    """
        delly call ${bam} -v ${merged_bcf} -g ${ref} -o ${isotype}_geno.bcf
        bcftools index -f ${isotype}_geno.bcf
    """

}

process proc_genos {

    label "dell_big"
        
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        path("*")

    output:
        tuple file("WI.DELLYpif.germline-filter.vcf.gz"), file("WI.DELLYpif.germline-filter.vcf.gz.tbi"), emit: delly_germline_filtered
        tuple file("WI.DELLYpif.raw.bed"), emit: raw_ceandr_bed
        tuple file("WI.DELLYpif.raw.vcf.gz"), file("WI.DELLYpif.raw.vcf.gz.tbi"), file("WI.DELLYpif.raw.stats.txt"), file("WI.DELLYpif.germline-filter.stats.txt")

    script:
        """
            ls *.bcf > bcf_list.txt

            bcftools merge -m id -Ob -o WI.DELLYpif.raw.bcf -l bcf_list.txt
            bcftools index -f WI.DELLYpif.raw.bcf

            bcftools query -l WI.DELLYpif.raw.bcf | sort > sample_names.txt

            bcftools view --samples-file=sample_names.txt -Oz -o WI.DELLYpif.raw.vcf.gz WI.DELLYpif.raw.bcf
            tabix -p vcf -f WI.DELLYpif.raw.vcf.gz

            delly filter -f germline WI.DELLYpif.raw.bcf -o WI.DELLYpif.germline-filter.bcf
            bcftools view -Oz -o WI.DELLYpif.germline-filter.vcf.gz WI.DELLYpif.germline-filter.bcf
            tabix -p vcf -f WI.DELLYpif.germline-filter.vcf.gz

            bcftools view --samples-file=sample_names.txt -Oz -o WI.DELLYpif.germline-filter.vcf.gz WI.DELLYpif.germline-filter.bcf
            tabix -p vcf -f WI.DELLYpif.germline-filter.vcf.gz

            #bcftools query WI.DELLYpif.germline-filter.vcf.gz -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' |\\
            #awk -F"|" '{print \$1, \$2, \$3, \$4, \$5}' OFS="\\t" > WI.DELLYpif.germline.bed
            
            ## I switched this process to send the raw SVs to the output_caendr_pif process, but this process still applies the delly filter command above
            ## THIS CAN BE MADE MORE EFFICIENT WHEN WE SETTLE ON WHICH FILE TO USE
            
            bcftools query WI.DELLYpif.raw.vcf.gz -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' |\\
            awk -F"|" '{print \$1, \$2, \$3, \$4, \$5}' OFS="\\t" > WI.DELLYpif.raw.bed
            
            bcftools stats --verbose WI.DELLYpif.raw.vcf.gz > WI.DELLYpif.raw.stats.txt
            bcftools stats --verbose WI.DELLYpif.germline-filter.vcf.gz > WI.DELLYpif.germline-filter.stats.txt     
        """

}

process proc_genos_all {

    label "dell_big"
        
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        path("*")

    output:
        tuple file("WI.DELLY.germline-filter.vcf.gz"), file("WI.DELLY.germline-filter.vcf.gz.tbi"), emit: delly_germline_filtered
      //  tuple file("WI.DELLY.raw.bed"), emit: raw_ceandr_bed
        tuple file("WI.DELLY.raw.vcf.gz"), file("WI.DELLY.raw.vcf.gz.tbi"), file("WI.DELLY.raw.stats.txt"), file("WI.DELLY.germline-filter.stats.txt")

    script:
        """
            ls *.bcf > bcf_list.txt

            bcftools merge -m id -Ob -o WI.DELLY.raw.bcf -l bcf_list.txt
            bcftools index -f WI.DELLY.raw.bcf

            bcftools query -l WI.DELLY.raw.bcf | sort > sample_names.txt

            bcftools view --samples-file=sample_names.txt -Oz -o WI.DELLY.raw.vcf.gz WI.DELLY.raw.bcf
            tabix -p vcf -f WI.DELLY.raw.vcf.gz

            delly filter -f germline WI.DELLY.raw.bcf -o WI.DELLY.germline-filter.bcf
            bcftools view -Oz -o WI.DELLY.germline-filter.vcf.gz WI.DELLY.germline-filter.bcf
            tabix -p vcf -f WI.DELLY.germline-filter.vcf.gz

            bcftools view --samples-file=sample_names.txt -Oz -o WI.DELLY.germline-filter.vcf.gz WI.DELLY.germline-filter.bcf
            tabix -p vcf -f WI.DELLY.germline-filter.vcf.gz

            #bcftools query WI.DELLY.germline-filter.vcf.gz -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' |\\
            #awk -F"|" '{print \$1, \$2, \$3, \$4, \$5}' OFS="\\t" > WI.DELLY.germline.bed
            
            ## I switched this process to send the raw SVs to the output_caendr_pif process, but this process still applies the delly filter command above
            ## THIS CAN BE MADE MORE EFFICIENT WHEN WE SETTLE ON WHICH FILE TO USE
            
            #bcftools query WI.DELLY.raw.vcf.gz -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\n]' |\\
            #awk -F"|" '{print \$1, \$2, \$3, \$4, \$5}' OFS="\\t" > WI.DELLY.raw.bed
            
            bcftools stats --verbose WI.DELLY.raw.vcf.gz > WI.DELLY.raw.stats.txt
            bcftools stats --verbose WI.DELLY.germline-filter.vcf.gz > WI.DELLY.germline-filter.stats.txt     
        """

}

process output_caendr_pif {

    publishDir "${params.out}/caendr_pif", mode: 'copy'

    label "R"
        
    input:
        tuple file(caendr_pif_script), file(bed), file(ref), val(species)

    output:
        tuple file("caendr.pif.${species}.bed.gz"), file("caendr.pif.${species}.bed.gz.tbi"), file("caendr.pif.${species}.vcf.gz"), file("caendr.pif.${species}.vcf.gz.csi"), emit: ceandr_pif

    """
        # Use config script to setup delly run parameters
        Rscript --vanilla ${caendr_pif_script} ${bed} ${ref} ${species}
    """

}