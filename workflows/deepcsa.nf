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
include { INPUT_CHECK                                       } from '../subworkflows/local/input_check'

include { DEPTH_ANALYSIS            as DEPTHANALYSIS        } from '../subworkflows/local/depthanalysis/main'
include { CREATE_PANELS             as CREATEPANELS         } from '../subworkflows/local/createpanels/main'

include { MUTATION_PREPROCESSING    as MUT_PREPROCESSING    } from '../subworkflows/local/mutationpreprocessing/main'

include { MUTATION_RATE          as MUTRATE           } from '../subworkflows/local/mutationrate/main'

include { MUTATIONAL_PROFILE        as MUTPROFILEALL        } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILENONPROT    } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILEEXONS      } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILEINTRONS    } from '../subworkflows/local/mutationprofile/main'

include { MUTABILITY                as MUTABILITYALL        } from '../subworkflows/local/mutability/main'
include { MUTABILITY                as MUTABILITYNONPROT    } from '../subworkflows/local/mutability/main'

include { ONCODRIVEFML_ANALYSIS     as ONCODRIVEFMLALL      } from '../subworkflows/local/oncodrivefml/main'
include { ONCODRIVEFML_ANALYSIS     as ONCODRIVEFMLNONPROT  } from '../subworkflows/local/oncodrivefml/main'
include { ONCODRIVE3D_ANALYSIS      as ONCODRIVE3D          } from '../subworkflows/local/oncodrive3d/main'
include { ONCODRIVECLUSTL_ANALYSIS  as ONCODRIVECLUSTL      } from '../subworkflows/local/oncodriveclustl/main'
// include { OMEGA_ANALYSIS            as OMEGA                } from '../subworkflows/local/omega/main'


include { SIGNATURES                as SIGNATURESALL        } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESNONPROT    } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESEXONS      } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESINTRONS    } from '../subworkflows/local/signatures/main'

// include { DEPTH_ANALYSIS            as MUTRATE              } from '../subworkflows/local/depthanalysis/main'

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

include { ANNOTATE_DEPTHS           as ANNOTATEDEPTHS           } from '../modules/local/annotatedepth/main'


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


    //
    // Separate input BAMs and VCFs
    //
    INPUT_CHECK.out.mutations.
    map{ it -> [it[0], it[1]]}.
    set{ meta_vcfs_alone }

    INPUT_CHECK.out.mutations.
    map{ it -> [it[0], it[2]]}.
    set{ meta_bams_alone }



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


    // Depth analysis: compute and plots
    DEPTHANALYSIS(meta_bams_alone)
    ch_versions = ch_versions.mix(DEPTHANALYSIS.out.versions)

    // Panels generation: all modalities
    CREATEPANELS(DEPTHANALYSIS.out.depths, vep_cache, vep_extra_files)
    ch_versions = ch_versions.mix(CREATEPANELS.out.versions)

    ANNOTATEDEPTHS(DEPTHANALYSIS.out.depths, CREATEPANELS.out.all_panel)
    ANNOTATEDEPTHS.out.annotated_depths.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ annotated_depths }

    // Mutation preprocessing
    MUT_PREPROCESSING(meta_vcfs_alone, vep_cache, vep_extra_files, CREATEPANELS.out.exons_consensus_bed)
    ch_versions = ch_versions.mix(MUT_PREPROCESSING.out.versions)


    // Mutation Rate
    // TODO: input all the bedfiles in the same channel
    consensus_annot_panel_all = Channel.of([ [ id: "consensus_panel_all" ], params.consensus_annot_panel_all ])
    consensus_annot_panel_protaffect = Channel.of([ [ id: "consensus_panel_protaffect" ], params.consensus_annot_panel_protaffect ])
    consensus_annot_panel_nonprotaffect = Channel.of([ [ id: "consensus_panel_nonprotaffect" ], params.consensus_annot_panel_nonprotaffect ])
    depths = Channel.of([ [ id: "all_samples" ], params.depths ])
    maf = Channel.of([ [ id: "all_samples" ], params.maf ])
    MUTRATE(maf, depths, consensus_annot_panel_all, consensus_annot_panel_protaffect, consensus_annot_panel_nonprotaffect)
    ch_versions = ch_versions.mix(MUTRATE.out.versions)

    // Mutational profile
    MUTPROFILEALL(MUT_PREPROCESSING.out.mafs, annotated_depths, CREATEPANELS.out.all_bed)
    MUTPROFILENONPROT(MUT_PREPROCESSING.out.mafs, annotated_depths, CREATEPANELS.out.nonprot_consensus_bed)
    MUTPROFILEEXONS(MUT_PREPROCESSING.out.mafs, annotated_depths, CREATEPANELS.out.exons_consensus_bed)
    MUTPROFILEINTRONS(MUT_PREPROCESSING.out.mafs, annotated_depths, CREATEPANELS.out.introns_consensus_bed)

    ch_versions = ch_versions.mix(MUTPROFILEALL.out.versions)

    MUTABILITYALL(MUT_PREPROCESSING.out.mafs,
                    annotated_depths,
                    MUTPROFILEALL.out.profile,
                    CREATEPANELS.out.exons_consensus_panel,
                    CREATEPANELS.out.exons_consensus_bed
                    )

    MUTABILITYNONPROT(MUT_PREPROCESSING.out.mafs,
                        annotated_depths,
                        MUTPROFILENONPROT.out.profile,
                        CREATEPANELS.out.exons_consensus_panel,
                        CREATEPANELS.out.exons_consensus_bed
                        )

    //
    // Positive selection
    //

    MUT_PREPROCESSING.out.mafs
    .join(MUTABILITYALL.out.mutability)
    .set{mutations_n_mutabilitiesall}

    MUT_PREPROCESSING.out.mafs
    .join(MUTABILITYNONPROT.out.mutability)
    .set{mutations_n_mutabilitiesnonprot}

    // OncodriveFML
    ONCODRIVEFMLALL(mutations_n_mutabilitiesall, CREATEPANELS.out.exons_consensus_panel)
    ONCODRIVEFMLNONPROT(mutations_n_mutabilitiesnonprot, CREATEPANELS.out.exons_consensus_panel)
    ch_versions = ch_versions.mix(ONCODRIVEFMLALL.out.versions)

    // Oncodrive3D
    ONCODRIVE3D(mutations_n_mutabilitiesall)
    ch_versions = ch_versions.mix(ONCODRIVE3D.out.versions)

    // Omega
    // MUT_PREPROCESSING.out.mafs
    // .join(MUTPROFILE.out.profile)
    // .set{mutations_n_profile}

    // mutations_n_profile
    // .join(depths)
    // .set{mutations_n_profile_n_depths}

    // annotated_panel should come from depth analysis
    // OMEGA(mutations_n_profile_n_depths, annotated_panel)

    // OncodriveClustl
    ONCODRIVECLUSTL(mutations_n_mutabilitiesall, CREATEPANELS.out.exons_consensus_panel)


    // Signature Analysis
    SIGNATURESALL(MUTPROFILEALL.out.wgs_sigprofiler, params.cosmic_ref_signatures)
    SIGNATURESNONPROT(MUTPROFILENONPROT.out.wgs_sigprofiler, params.cosmic_ref_signatures)
    SIGNATURESEXONS(MUTPROFILEEXONS.out.wgs_sigprofiler, params.cosmic_ref_signatures)
    SIGNATURESINTRONS(MUTPROFILEINTRONS.out.wgs_sigprofiler, params.cosmic_ref_signatures)
    ch_versions = ch_versions.mix(SIGNATURESALL.out.versions)



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
