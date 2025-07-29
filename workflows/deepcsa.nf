/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_deepcsa'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

include { PLOT_DEPTHS               as PLOTDEPTHSALLCONS    } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONS      } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONSCONS  } from '../subworkflows/local/plotdepths/main'

include { MUTATION_PREPROCESSING    as MUT_PREPROCESSING    } from '../subworkflows/local/mutationpreprocessing/main'

include { MUTATION_DENSITY          as MUTDENSITYALL           } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYPROT          } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYNONPROT       } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYSYNONYMOUS    } from '../subworkflows/local/mutationdensity/main'

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

include { INDELS_SELECTION          as INDELSSELECTION      } from '../subworkflows/local/indels/main'

include { MUTATED_CELLS_VAF         as MUTATEDCELLSVAF      } from '../subworkflows/local/mutatedcells/vaf/main'

include { EXPECTED_MUTATED_CELLS    as EXPECTEDMUTATEDCELLS } from '../subworkflows/local/mutatedcells/expected/main'

include { SIGNATURES                as SIGNATURESALL        } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESNONPROT    } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESEXONS      } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESINTRONS    } from '../subworkflows/local/signatures/main'

include { PLOTTING_SUMMARY          as PLOTTINGSUMMARY      } from '../subworkflows/local/plottingsummary/main'

include { REGRESSIONS               as REGRESSIONSMUTDENSITY          } from '../subworkflows/local/regressions/main'
include { REGRESSIONS               as REGRESSIONSONCODRIVEFML     } from '../subworkflows/local/regressions/main'
include { REGRESSIONS               as REGRESSIONSOMEGA            } from '../subworkflows/local/regressions/main'
include { REGRESSIONS               as REGRESSIONSOMEGAGLOB        } from '../subworkflows/local/regressions/main'

include { DNDS                      as DNDS                 } from '../subworkflows/local/dnds/main'

include { TABIX_BGZIPTABIX_QUERY    as DEPTHSALLCONS        } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSEXONSCONS      } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSPROTCONS       } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSNONPROTCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSINTRONSCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSSYNONYMOUSCONS } from '../modules/nf-core/tabix/bgziptabixquery/main'

include { SELECT_MUTDENSITIES       as SYNMUTDENSITY        } from '../modules/local/select_mutdensity/main'
include { SELECT_MUTDENSITIES       as SYNMUTREADSRATE      } from '../modules/local/select_mutdensity/main'

include { DNA_2_PROTEIN_MAPPING     as DNA2PROTEINMAPPING   } from '../modules/local/dna2protein/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installed directly from nf-core/modules.
*/

// MODULE
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { ANNOTATE_DEPTHS           as ANNOTATEDEPTHS       } from '../modules/local/annotatedepth/main'
include { DOWNSAMPLE_DEPTHS         as DOWNSAMPLEDEPTHS     } from '../modules/local/downsample/depths/main'
include { TABLE_2_GROUP             as TABLE2GROUP          } from '../modules/local/table2groups/main'
include { MUTATIONS_2_SIGNATURES    as MUTS2SIGS            } from '../modules/local/mutations2sbs/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEEPCSA{

    // // Input channel definitions
    features_table  = Channel.fromPath( params.features_table ?: params.input, checkIfExists: true)

    wgs_trinucs     = params.wgs_trinuc_counts
                            ? Channel.fromPath( params.wgs_trinuc_counts, checkIfExists: true).first()
                            : Channel.empty()
    cosmic_ref      = params.cosmic_ref_signatures
                            ? Channel.fromPath( params.cosmic_ref_signatures, checkIfExists: true).first()
                            : Channel.empty()
    datasets3d      = params.datasets3d
                            ? Channel.fromPath( params.datasets3d, checkIfExists: true).first()
                            : Channel.empty()
    annotations3d   = params.annotations3d
                            ? Channel.fromPath( params.annotations3d, checkIfExists: true).first()
                            : Channel.empty()
    seqinfo_df      = params.datasets3d
                            ? Channel.fromPath( "${params.datasets3d}/seq_for_mut_prob.tsv", checkIfExists: true).first()
                            : Channel.empty()
    cadd_scores     = params.cadd_scores
                            ? Channel.of([
                                file(params.cadd_scores, checkIfExists : true),
                                file(params.cadd_scores_ind, checkIfExists : true)
                                ]).first()
                            : Channel.empty()


    // if the user wants to use custom gene groups, import the gene groups table
    // otherwise I am using the input csv as a dummy value channel
    custom_groups_table = params.custom_groups_file
                                ? Channel.fromPath( params.custom_groups_file, checkIfExists: true).first()
                                : Channel.fromPath(params.input)

    // if the user wants to use custom BED file for computing the depths, import the BED file
    // otherwise I am using the input csv as a dummy value channel
    custom_bed_file     = params.custom_bedfile
                                ? Channel.fromPath( params.custom_bedfile, checkIfExists: true).first()
                                : Channel.fromPath(params.input)

    // Initialize booleans based on user params
    def run_mutabilities    = (params.oncodrivefml || params.oncodriveclustl || params.oncodrive3d)
    def run_mutdensity         = (params.mutationdensity || params.omega)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )


    //
    // Separate input BAMs and VCFs
    //
    INPUT_CHECK.out.mutations
    .map{ it -> [ "id" : it[0].id ]}
    .set{ meta_samples_alone }

    INPUT_CHECK.out.mutations
    .map{ it -> [it[0], it[1]]}
    .set{ meta_vcfs_alone }

    INPUT_CHECK.out.mutations
    .map{ it -> [it[0], it[2]]}
    .set{ meta_bams_alone }


    TABLE2GROUP(features_table)


    // Depths and panel creation should be a single subworkflow
    // Depth analysis: compute and plots
    DEPTHANALYSIS(meta_bams_alone, custom_bed_file)

    // Panels generation: all modalities
    CREATEPANELS(DEPTHANALYSIS.out.depths)

    ANNOTATEDEPTHS(DEPTHANALYSIS.out.depths, CREATEPANELS.out.all_panel, TABLE2GROUP.out.json_allgroups, file(params.input))
    ANNOTATEDEPTHS.out.annotated_depths.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ annotated_depths_full }

    // if (params.downsample && params.downsample_proportion < 1) {
    if (params.downsample ){
        DOWNSAMPLEDEPTHS(annotated_depths_full)
        annotated_depths = DOWNSAMPLEDEPTHS.out.downsampled_depths
    } else {
        annotated_depths = annotated_depths_full
    }

    PLOTDEPTHSALLCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
    PLOTDEPTHSEXONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_bed, CREATEPANELS.out.exons_panel)
    PLOTDEPTHSEXONSCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_consensus_bed, CREATEPANELS.out.exons_consensus_panel)

    // Mutation preprocessing
    MUT_PREPROCESSING(meta_vcfs_alone,
                        CREATEPANELS.out.all_consensus_bed,
                        CREATEPANELS.out.exons_bed,
                        TABLE2GROUP.out.json_allgroups,
                        seqinfo_df,
                        CREATEPANELS.out.added_custom_regions
                        )
    somatic_mutations = MUT_PREPROCESSING.out.somatic_mafs

    positive_selection_results = somatic_mutations


    if (params.vep_species == 'homo_sapiens'){
        DNA2PROTEINMAPPING(MUT_PREPROCESSING.out.mutations_all_samples, CREATEPANELS.out.exons_consensus_panel)
    }


    // Intersect BED of desired sites with samples' depths
    DEPTHSALLCONS(annotated_depths, CREATEPANELS.out.all_consensus_bed)
    DEPTHSEXONSCONS(annotated_depths, CREATEPANELS.out.exons_consensus_bed)
    if (run_mutdensity || params.profilenonprot) {
        DEPTHSNONPROTCONS(annotated_depths, CREATEPANELS.out.nonprot_consensus_bed)
    }

    if (run_mutdensity) {
        DEPTHSPROTCONS(annotated_depths, CREATEPANELS.out.prot_consensus_bed)
        DEPTHSSYNONYMOUSCONS(annotated_depths, CREATEPANELS.out.synonymous_consensus_bed)
    }


    if (run_mutdensity){
        // Mutation Density
        MUTDENSITYALL(somatic_mutations, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
        MUTDENSITYPROT(somatic_mutations, DEPTHSPROTCONS.out.subset, CREATEPANELS.out.prot_consensus_bed, CREATEPANELS.out.prot_consensus_panel)
        MUTDENSITYNONPROT(somatic_mutations, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, CREATEPANELS.out.nonprot_consensus_panel)
        MUTDENSITYSYNONYMOUS(somatic_mutations, DEPTHSSYNONYMOUSCONS.out.subset, CREATEPANELS.out.synonymous_consensus_bed, CREATEPANELS.out.synonymous_consensus_panel)

        Channel.of([ [ id: "all_samples" ] ])
        .join( MUTDENSITYSYNONYMOUS.out.mutdensities )
        .set{ all_samples_syn_mutdensity }

        SYNMUTDENSITY(all_samples_syn_mutdensity)

        SYNMUTREADSRATE(all_samples_syn_mutdensity)


        // Concatenate all outputs into a single file
        Channel.empty()
        .concat(MUTDENSITYALL.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYPROT.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYNONPROT.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYSYNONYMOUS.out.mutdensities.map{ it -> it[1]}.flatten())
        .set{ all_mutdensities }
        all_mutdensities.collectFile(name: "all_mutdensities.tsv", storeDir:"${params.outdir}/mutdensity", skip: 1, keepHeader: true).set{ all_mutdensities_file }

    }


    // Mutational profile
    if (params.profileall){
        MUTPROFILEALL(somatic_mutations, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, wgs_trinucs)
    }
    if (params.profilenonprot){
        MUTPROFILENONPROT(somatic_mutations, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, wgs_trinucs)
    }
    if (params.profileexons){
        MUTPROFILEEXONS(somatic_mutations, DEPTHSEXONSCONS.out.subset, CREATEPANELS.out.exons_consensus_bed, wgs_trinucs)
    }
    if (params.profileintrons){
        DEPTHSINTRONSCONS(annotated_depths, CREATEPANELS.out.introns_consensus_bed)
        MUTPROFILEINTRONS(somatic_mutations, DEPTHSINTRONSCONS.out.subset, CREATEPANELS.out.introns_consensus_bed, wgs_trinucs)
    }


    if (run_mutabilities) {
        if (params.profileall){
            MUTABILITYALL(somatic_mutations,
                            annotated_depths,
                            MUTPROFILEALL.out.profile,
                            CREATEPANELS.out.exons_consensus_panel,
                            CREATEPANELS.out.exons_consensus_bed
                            )
        }
        if (params.profilenonprot){
            MUTABILITYNONPROT(somatic_mutations,
                                annotated_depths,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_panel,
                                CREATEPANELS.out.exons_consensus_bed
                                )
        }
    }


    if (params.indels){
        INDELSSELECTION(somatic_mutations,
                        CREATEPANELS.out.all_consensus_bed
                        )
        positive_selection_results = positive_selection_results.join(INDELSSELECTION.out.indels, remainder: true)
    }

    if (params.expected_mutated_cells){
        EXPECTEDMUTATEDCELLS(MUT_PREPROCESSING.out.mutations_all_samples,
                                CREATEPANELS.out.exons_consensus_bed,
                                CREATEPANELS.out.exons_consensus_panel,
                                ANNOTATEDEPTHS.out.all_samples_depths,
                                CREATEPANELS.out.full_panel_annotated
                                )
    }

    //
    // Positive selection
    //

    // OncodriveFML
    if (params.oncodrivefml){
        oncodrivefml_regressions_files = Channel.empty()
        if (params.profileall){
            mode = "all"
            ONCODRIVEFMLALL(somatic_mutations, MUTABILITYALL.out.mutability,
                                CREATEPANELS.out.exons_consensus_panel,
                                cadd_scores, mode
                            )
            positive_selection_results = positive_selection_results.join(ONCODRIVEFMLALL.out.results_snvs, remainder: true)

            if (params.regressions){
                oncodrivefml_regressions_files = oncodrivefml_regressions_files.mix(ONCODRIVEFMLALL.out.results_snvs_folder.map{ it -> it[1] })
            }
        }
        if (params.profilenonprot){
            mode = "non_prot_aff"
            ONCODRIVEFMLNONPROT(somatic_mutations, MUTABILITYNONPROT.out.mutability,
                                    CREATEPANELS.out.exons_consensus_panel,
                                    cadd_scores, mode
                                )
            if (params.regressions){
                oncodrivefml_regressions_files = oncodrivefml_regressions_files.mix(ONCODRIVEFMLNONPROT.out.results_snvs_folder.map{ it -> it[1] })
            }
        }
    }

    if (params.oncodrive3d){
        if (params.profileall){
            // Oncodrive3D
            ONCODRIVE3D(somatic_mutations, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_bed,
                        datasets3d, annotations3d, MUT_PREPROCESSING.out.all_raw_vep_annotation)
            positive_selection_results = positive_selection_results.join(ONCODRIVE3D.out.results, remainder: true)
            positive_selection_results = positive_selection_results.join(ONCODRIVE3D.out.results_pos, remainder: true)

        }
    }

    // if (params.expected_mutated_cells & params.dnds){
    if (params.dnds){
        covariates = params.dnds_covariates ? Channel.fromPath( params.dnds_covariates, checkIfExists: true).first() : Channel.empty()
        ref_transcripts = params.dnds_ref_transcripts ? Channel.fromPath( params.dnds_ref_transcripts, checkIfExists: true).first() : Channel.empty()
        DNDS(somatic_mutations,
                    DEPTHSEXONSCONS.out.subset,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel,
                    covariates,
                    ref_transcripts
                    )
    }

    if (params.omega){
        omega_regressions_files = Channel.empty()
        omega_regressions_files_gloc = Channel.empty()

        // Omega
        if (params.profileall){
            OMEGA(somatic_mutations,
                    DEPTHSEXONSCONS.out.subset,
                    MUTPROFILEALL.out.profile,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel,
                    custom_groups_table,
                    CREATEPANELS.out.domains_panel_bed,
                    SYNMUTDENSITY.out.mutdensity,
                    CREATEPANELS.out.panel_annotated_rich
                    )
            positive_selection_results = positive_selection_results.join(OMEGA.out.results, remainder: true)
            positive_selection_results = positive_selection_results.join(OMEGA.out.results_global, remainder: true)
            site_comparison_results = OMEGA.out.site_comparison

            if (params.regressions){
                omega_regressions_files = omega_regressions_files.mix(OMEGA.out.results.map{ it -> it[1] })
                omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGA.out.results_global.map{ it -> it[1] })
            }

            if (params.omega_multi){
                // Omega multi
                OMEGAMULTI(somatic_mutations,
                            DEPTHSEXONSCONS.out.subset,
                            MUTPROFILEALL.out.profile,
                            CREATEPANELS.out.exons_consensus_bed,
                            CREATEPANELS.out.exons_consensus_panel,
                            custom_groups_table,
                            CREATEPANELS.out.domains_panel_bed,
                            SYNMUTREADSRATE.out.mutdensity,
                            CREATEPANELS.out.panel_annotated_rich
                            )
                positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results, remainder: true)
                positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results_global, remainder: true)
                if (params.regressions){
                    omega_regressions_files = omega_regressions_files.mix(OMEGAMULTI.out.results.map{ it -> it[1] })
                    omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGAMULTI.out.results_global.map{ it -> it[1] })
                }
            }
        }
        if (params.profilenonprot){
            OMEGANONPROT(somatic_mutations,
                            DEPTHSEXONSCONS.out.subset,
                            MUTPROFILENONPROT.out.profile,
                            CREATEPANELS.out.exons_consensus_bed,
                            CREATEPANELS.out.exons_consensus_panel,
                            custom_groups_table,
                            CREATEPANELS.out.domains_panel_bed,
                            SYNMUTDENSITY.out.mutdensity,
                            CREATEPANELS.out.panel_annotated_rich
                            )
            if (params.regressions){
                omega_regressions_files = omega_regressions_files.mix(OMEGANONPROT.out.results.map{ it -> it[1] })
                omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGANONPROT.out.results_global.map{ it -> it[1] })
            }

            if (params.omega_multi){
                OMEGANONPROTMULTI(somatic_mutations,
                                    DEPTHSEXONSCONS.out.subset,
                                    MUTPROFILENONPROT.out.profile,
                                    CREATEPANELS.out.exons_consensus_bed,
                                    CREATEPANELS.out.exons_consensus_panel,
                                    custom_groups_table,
                                    CREATEPANELS.out.domains_panel_bed,
                                    SYNMUTREADSRATE.out.mutdensity,
                                    CREATEPANELS.out.panel_annotated_rich
                                    )

                if (params.regressions){
                    omega_regressions_files = omega_regressions_files.mix(OMEGANONPROTMULTI.out.results.map{ it -> it[1] })
                    omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGANONPROTMULTI.out.results_global.map{ it -> it[1] })
                }
            }

        }

    }

    if (params.mutated_cells_vaf){
        MUT_PREPROCESSING.out.somatic_mafs
        .join(meta_samples_alone)
        .set{ sample_mutations_only }

        MUTATEDCELLSVAF(sample_mutations_only,
                                CREATEPANELS.out.exons_consensus_bed,
                                OMEGA.out.results_global,
                                features_table
                                // OMEGAMULTI.out.results_global
                                )
    }


    if (params.oncodriveclustl){
        // OncodriveClustl
        if (params.profileall){
            ONCODRIVECLUSTL(somatic_mutations, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_panel)
        }
    }


    if (params.signatures){

        // Signature Analysis
        if (params.profileall){
            SIGNATURESALL(MUTPROFILEALL.out.wgs_sigprofiler, cosmic_ref, TABLE2GROUP.out.json_samples)

            somatic_mutations
            .join(SIGNATURESALL.out.mutation_probs)
            .set{mutations_n_sbs}
            MUTS2SIGS(mutations_n_sbs)


        }
        if (params.profilenonprot){
            SIGNATURESNONPROT(MUTPROFILENONPROT.out.wgs_sigprofiler, cosmic_ref, TABLE2GROUP.out.json_samples)
        }
        if (params.profileexons){
            SIGNATURESEXONS(MUTPROFILEEXONS.out.wgs_sigprofiler, cosmic_ref, TABLE2GROUP.out.json_samples)
        }
        if (params.profileintrons){
            SIGNATURESINTRONS(MUTPROFILEINTRONS.out.wgs_sigprofiler, cosmic_ref, TABLE2GROUP.out.json_samples)
        }
    }

    if ( params.omega & params.oncodrive3d & params.oncodrivefml & params.indels ){
        positive_selection_results_ready = positive_selection_results.map { element -> [element[0], element[1..-1]] }
        PLOTTINGSUMMARY(positive_selection_results_ready,
                        all_mutdensities_file,
                        site_comparison_results,
                        ANNOTATEDEPTHS.out.all_samples_depths,
                        TABLE2GROUP.out.json_samples,
                        TABLE2GROUP.out.json_allgroups,
                        CREATEPANELS.out.exons_consensus_panel,
                        CREATEPANELS.out.panel_annotated_rich,
                        seqinfo_df
                        )
    }

    // Regressions
    if (params.regressions){

        if (params.mutationdensity && params.mutdensity_regressions){
            REGRESSIONSMUTDENSITY("mutdensity", all_mutdensities_file, params.mutdensity_regressions)
        }

        if (params.oncodrivefml && params.oncodrivefml_regressions){
            REGRESSIONSONCODRIVEFML("oncodrivefml", oncodrivefml_regressions_files.toList(), params.oncodrivefml_regressions)
        }

        if (params.omega && params.omega_regressions){
            REGRESSIONSOMEGA("omega", omega_regressions_files.toList(), params.omega_regressions)
        }

        if (params.omega_globalloc && params.omega_regressions){
            REGRESSIONSOMEGAGLOB("omegagloballoc", omega_regressions_files_gloc.toList(), params.omega_regressions)
        }

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        Channel.topic('versions').unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // Info required for completion email and summary
    def multiqc_report = []

    def summary_params = paramsSummaryMap(workflow)
    workflow_summary    = paramsSummaryMultiqc(summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = methodsDescriptionText(ch_multiqc_custom_methods_description)
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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
