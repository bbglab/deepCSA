/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDeepcsa.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Consisting of a mix of local and nf-core/modules.
*/

// SUBWORKFLOW
include { INPUT_CHECK                                 } from '../subworkflows/local/input_check'
include { ONCODRIVEFML_ANALYSIS  as ONCODRIVEFML      } from '../subworkflows/local/oncodrivefml/main'
include { ONCODRIVE3D_ANALYSIS   as ONCODRIVE3D       } from '../subworkflows/local/oncodrive3d/main'
include { MUTATION_PREPROCESSING as MUT_PREPROCESSING } from '../subworkflows/local/mutationpreprocessing/main'

// include { DEPTH_ANALYSIS as DEPTHANALYSIS       } from '../subworkflows/local/depthanalysis/main'

// Download annotation cache if needed
include { PREPARE_CACHE                               } from '../subworkflows/local/prepare_cache/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installed directly from nf-core/modules.
*/

// MODULE
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DEEPCSA {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema


    INPUT_CHECK.out.mutations.
    map{ it -> [it[0], it[1]]}.
    set{ meta_vcfs_alone }

    INPUT_CHECK.out.mutations.
    map{ it -> [it[0], it[2]]}.
    set{ meta_bams_alone }

    // // Run depth analysis subworkflow
    // DEPTHANALYSIS(INPUT_CHECK.mutations)


    // TODO: test if downloading VEP cache works
    // Download Ensembl VEP cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache ? [] : Channel.of([ [ id:"${params.vep_genome}.${params.vep_cache_version}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info)
        vep_cache = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
        ch_versions = ch_versions.mix(PREPARE_CACHE.out.versions)
    } else {
        vep_cache = params.vep_cache
    }
    vep_extra_files = []


    bedfile = params.bedf
    MUT_PREPROCESSING(meta_vcfs_alone, vep_cache, vep_extra_files, bedfile)
    ch_versions = ch_versions.mix(MUT_PREPROCESSING.out.versions)

    // ONCODRIVEFML(params.muts, params.mutabs, params.mutabs_index, params.bedf)

    // ONCODRIVE3D(params.muts_3d, params.mutabs, params.mutabs_index)


// SUBWORKFLOWS
//     depths analysis
//     mutation preprocessing
//     mutational profile
//     Module to subsets a MAF

//     oncodrivefml
//     omega
//     oncodrive3d
//     oncodriveclustl

//     mutation rate

//     signature extraction?
//     signature decomposition of our profiles






// Main containers:
//     - pybedtools                        avail
//     - pandas and basic plotting         avail
//     - bgreference                                   ¡¡highest priority!!
//     - SigProfilerMatrixGenerator        avail
//     - sigprofilerextractor
//     - omega
//     - oncodriveclustl                   ?local
//     - oncodrive3d                       local
//     - SigLasso





//     // SUBWORKFLOW: depths analysis
//     //     Define the regions to analyse (extended regions)
//     //         Input:
//     //             Target BED file?
//     //             Depths files, load all the depths files per sample.

//     //     Read the input depth files -> Build depths matrix, index with tabix maybe?

//     //         Output:

//     //         BED file

//     //         Report information of the regions in terms of:
//     //             % exons
//     //             % introns
//     //             % intergenic / off target (if any)
//     //                 any specific region with high coverage

//     //             Undercovered exons

//     //             Overall depth tendency per region

//     //             Heterogeneity between samples per region

//     // Modules list:
//         - Build depth matrix
//         - Estimate extended regions
//         - Define poorly covered regions -> provide a BED file of potential blacklists.
//         - ?Find depth correlated samples
//         - ?Find well covered off target regions






// Module that subsets a MAF to the format of interest in terms of:
//     - VAF
//     - Mutations passing certain filters
//     - Columns
//     - File format for: OncodriveFML, omega, oncodrive3d ...

// This can be used at the level of workflow or inside a subworkflow



// SUBWORKFLOW: mutational profile
// Description:
//     - Compute mutational profile in the regions of interest.
// Input:
//     - Depths
//     - Mutations
//     - Regions where mutations "can occur" (we could potentially run this in three different ways:
//                                                 with all sites
//                                                     protein-affecting sites,
//                                                     non-protein affecting
//                                             )
// Output:
//     - Tables with samples as columns and rows as contexts
//         96 mutational profile       ***
//         384 mutational profile
//         1536 mutational profile
//     - SigProfilerMatrixGenerator output?
//         - The advantage is that they compute the potential transcription strand bias. (384)
//     - Plots

// Modules:
//     - Intersect regions bedfile with depths
//     - Intersect the intersection of the previous two files with the mutations
//     - Compute mutational profiles per sample
//     - Compute matrix of mutations per context
//     - Extrapolate mutational profiles in the regions under analysis to whole genome trinucleotide or pentanucelotide contexts.







// SUBWORKFLOW: oncodrivefml
// Description:
//     - Compute positive selection based on functional impact bias.
// Input:
//     - Depths
//     - Mutations
//     - Mutation profile
//     - Regions to focus the analysis (we could potentially run this in different ways:
//                                         with all genes independently
//                                         grouping certain genes (this is not a priority at all!!!)
//                                         )
// Output:
//     - OncodriveFML results
//     - Plots ?

// Modules:
//     - Intersect regions bedfile with depths
//     - Intersect the intersection of the previous two files with the mutations
//     - Preprocessing module to compute mutabilities from mutational profile and depths.
//     - ?Preprocessing module to compute mutabilities from mutational profile and depths, and observed mutations per gene. (this is not a priority at all!!!)
//     - Run OncodriveFML for each sample individually.
//                                 (we could potentially run this for:
//                                         all samples
//                                         specific groups of samples (see provided to omega)
//                                         )
//     - Plot results:
//             Ideas:
//                 Clustermap of all individual samples




// SUBWORKFLOW: omega
// Description:
//     - Compute positive selection based on mutation recurrence.
// Input:
//     - Depths
//     - Mutations
//     - Mutation profile
//     - Regions to focus the analysis annotated (pending to update omega version with this working)
//     - Sample gene and consequence groups. (set defaults: all samples independent, all samples together
//                                                             all genes independent, all genes together,
//                                                                         all genes in the analysis that are known cancer genes,
//                                                             set default consequence types
//                                                             )
// Output:
//     - omega results
//     - Plots ?

// Modules:
//     - Intersect regions bedfile with depths
//     - Intersect the intersection of the previous two files with the mutations
//     - Omega preprocessing module
//     - Omega estimator module
//     - Plot results:
//             Ideas:
//                 Clustermap of all individual samples
//                 Sanger-like, with number of mutations per type with all samples, and then positive selection for each of the impacts
//                 Positive selection per each of the genes.


// Subworkflow signature assignment (or extraction)

//     - SigProfilerAssingment
//     - SigLasso
//     - MSA ( https://gitlab.com/s.senkin/MSA/-/tree/master )


//     - SigProfilerExtractor











    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDeepcsa.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDeepcsa.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
