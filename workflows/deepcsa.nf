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

features_table = params.features_table ? Channel.fromPath( params.features_table, checkIfExists: true) : Channel.fromPath(params.input)

wgs_trinucs = params.wgs_trinuc_counts ? Channel.fromPath( params.wgs_trinuc_counts, checkIfExists: true).first() : Channel.empty()
cosmic_ref = params.cosmic_ref_signatures ? Channel.fromPath( params.cosmic_ref_signatures, checkIfExists: true).first() : Channel.empty()
datasets3d = params.datasets3d ? Channel.fromPath( params.datasets3d, checkIfExists: true).first() : Channel.empty()


def run_mutabilities = (params.oncodrivefml || params.oncodriveclustl || params.oncodrive3d)


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

include { PLOT_DEPTHS               as PLOTDEPTHS           } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONS      } from '../subworkflows/local/plotdepths/main'

include { MUTATION_PREPROCESSING    as MUT_PREPROCESSING    } from '../subworkflows/local/mutationpreprocessing/main'

include { MUTATION_RATE             as MUTRATEALL           } from '../subworkflows/local/mutationrate/main'
include { MUTATION_RATE             as MUTRATEPROT          } from '../subworkflows/local/mutationrate/main'
include { MUTATION_RATE             as MUTRATENONPROT       } from '../subworkflows/local/mutationrate/main'

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

include { OMEGA_ANALYSIS            as OMEGA                } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGANONPROT         } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGAMULTI           } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGANONPROTMULTI    } from '../subworkflows/local/omega/main'



include { SIGNATURES                as SIGNATURESALL        } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESNONPROT    } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESEXONS      } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESINTRONS    } from '../subworkflows/local/signatures/main'

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
include { TABLE_2_GROUP             as TABLE2GROUP              } from '../modules/local/table2groups/main'
include { MUTATIONS_2_SIGNATURES    as MUTS2SIGS                } from '../modules/local/mutations2sbs/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DEEPCSA{

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


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


    TABLE2GROUP(features_table)
    ch_versions = ch_versions.mix(TABLE2GROUP.out.versions)


    // Depth analysis: compute and plots
    DEPTHANALYSIS(meta_bams_alone)
    ch_versions = ch_versions.mix(DEPTHANALYSIS.out.versions)

    // Panels generation: all modalities
    CREATEPANELS(DEPTHANALYSIS.out.depths, vep_cache, vep_extra_files)
    ch_versions = ch_versions.mix(CREATEPANELS.out.versions)

    ANNOTATEDEPTHS(DEPTHANALYSIS.out.depths, CREATEPANELS.out.all_panel, TABLE2GROUP.out.json_allgroups)
    ch_versions = ch_versions.mix(ANNOTATEDEPTHS.out.versions)
    ANNOTATEDEPTHS.out.annotated_depths.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ annotated_depths }

    PLOTDEPTHS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
    PLOTDEPTHSEXONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_bed, CREATEPANELS.out.exons_panel)

    // Mutation preprocessing
    MUT_PREPROCESSING(meta_vcfs_alone, vep_cache, vep_extra_files, CREATEPANELS.out.exons_consensus_bed, TABLE2GROUP.out.json_allgroups)
    ch_versions = ch_versions.mix(MUT_PREPROCESSING.out.versions)


    if (params.mutationrate){
        // Mutation Rate
        MUTRATEALL(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
        MUTRATEPROT(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.prot_consensus_bed, CREATEPANELS.out.prot_consensus_panel)
        MUTRATENONPROT(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.nonprot_consensus_bed, CREATEPANELS.out.nonprot_consensus_panel)
        ch_versions = ch_versions.mix(MUTRATEALL.out.versions)
        ch_versions = ch_versions.mix(MUTRATEPROT.out.versions)
        ch_versions = ch_versions.mix(MUTRATENONPROT.out.versions)

        // Concatenate all outputs into a single file
        mutrate_empty = Channel.empty()
        mutrate_empty
        .concat(MUTRATEALL.out.mutrates.map{ it -> it[1]}.flatten())
        .concat(MUTRATEPROT.out.mutrates.map{ it -> it[1]}.flatten())
        .concat(MUTRATENONPROT.out.mutrates.map{ it -> it[1]}.flatten())
        .set{ all_mutrates }
        all_mutrates.collectFile(name: "all_mutrates.tsv", storeDir:"${params.outdir}/mutrate", skip: 1, keepHeader: true)

    }


    // Mutational profile
    if (params.profileall){
        MUTPROFILEALL(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.all_bed, wgs_trinucs)
        ch_versions = ch_versions.mix(MUTPROFILEALL.out.versions)
    }
    if (params.profilenonprot){
        MUTPROFILENONPROT(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.nonprot_consensus_bed, wgs_trinucs)
        ch_versions = ch_versions.mix(MUTPROFILENONPROT.out.versions)
    }
    if (params.profileexons){
        MUTPROFILEEXONS(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.exons_consensus_bed, wgs_trinucs)
        ch_versions = ch_versions.mix(MUTPROFILEEXONS.out.versions)
    }
    if (params.profileintrons){
        MUTPROFILEINTRONS(MUT_PREPROCESSING.out.somatic_mafs, annotated_depths, CREATEPANELS.out.introns_consensus_bed, wgs_trinucs)
        ch_versions = ch_versions.mix(MUTPROFILEINTRONS.out.versions)
    }


    if (run_mutabilities) {
        if (params.profileall){
            MUTABILITYALL(MUT_PREPROCESSING.out.somatic_mafs,
                            annotated_depths,
                            MUTPROFILEALL.out.profile,
                            CREATEPANELS.out.exons_consensus_panel,
                            CREATEPANELS.out.exons_consensus_bed
                            )
            ch_versions = ch_versions.mix(MUTABILITYALL.out.versions)
        }
        if (params.profilenonprot){
            MUTABILITYNONPROT(MUT_PREPROCESSING.out.somatic_mafs,
                                annotated_depths,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_panel,
                                CREATEPANELS.out.exons_consensus_bed
                                )
            ch_versions = ch_versions.mix(MUTABILITYNONPROT.out.versions)
        }
    }



    //
    // Positive selection
    //

    // OncodriveFML
    if (params.oncodrivefml){
        if (params.profileall){
            ONCODRIVEFMLALL(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(ONCODRIVEFMLALL.out.versions)
        }
        if (params.profilenonprot){
            ONCODRIVEFMLNONPROT(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYNONPROT.out.mutability, CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(ONCODRIVEFMLNONPROT.out.versions)
        }
    }

    if (params.oncodrive3d){
        if (params.profileall){
            // Oncodrive3D
            ONCODRIVE3D(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_bed, datasets3d)
            ch_versions = ch_versions.mix(ONCODRIVE3D.out.versions)
        }
    }



    if (params.omega){

        // Omega
        if (params.profileall){
            OMEGA(MUT_PREPROCESSING.out.somatic_mafs,
                    annotated_depths,
                    MUTPROFILEALL.out.profile,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(OMEGA.out.versions)

            // Omega multi
            OMEGAMULTI(MUT_PREPROCESSING.out.somatic_mafs,
                        annotated_depths,
                        MUTPROFILEALL.out.profile,
                        CREATEPANELS.out.exons_consensus_bed,
                        CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(OMEGAMULTI.out.versions)
        }
        if (params.profilenonprot){
            OMEGANONPROT(MUT_PREPROCESSING.out.somatic_mafs,
                            annotated_depths,
                            MUTPROFILENONPROT.out.profile,
                            CREATEPANELS.out.exons_consensus_bed,
                            CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(OMEGANONPROT.out.versions)

            OMEGANONPROTMULTI(MUT_PREPROCESSING.out.somatic_mafs,
                                annotated_depths,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_bed,
                                CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(OMEGANONPROTMULTI.out.versions)
        }

    }

    if (params.oncodriveclustl){
        // OncodriveClustl
        if (params.profileall){
            ONCODRIVECLUSTL(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_panel)
            ch_versions = ch_versions.mix(ONCODRIVECLUSTL.out.versions)
        }
    }


    if (params.signatures){

        // Signature Analysis
        if (params.profileall){
            SIGNATURESALL(MUTPROFILEALL.out.matrix_sigprof, MUTPROFILEALL.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
            ch_versions = ch_versions.mix(SIGNATURESALL.out.versions)

            MUT_PREPROCESSING.out.somatic_mafs
            .join(SIGNATURESALL.out.mutation_probs)
            .set{mutations_n_sbs}
            MUTS2SIGS(mutations_n_sbs)
            ch_versions = ch_versions.mix(MUTS2SIGS.out.versions)


        }
        if (params.profilenonprot){
            SIGNATURESNONPROT(MUTPROFILENONPROT.out.matrix_sigprof, MUTPROFILENONPROT.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
            ch_versions = ch_versions.mix(SIGNATURESNONPROT.out.versions)
        }
        if (params.profileexons){
            SIGNATURESEXONS(MUTPROFILEEXONS.out.matrix_sigprof, MUTPROFILEEXONS.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
            ch_versions = ch_versions.mix(SIGNATURESEXONS.out.versions)
        }
        if (params.profileintrons){
            SIGNATURESINTRONS(MUTPROFILEINTRONS.out.matrix_sigprof, MUTPROFILEINTRONS.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
            ch_versions = ch_versions.mix(SIGNATURESINTRONS.out.versions)
        }


    }



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
