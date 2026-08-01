

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_deepcsa'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOW
include { INPUT_CHECK                                       } from '../subworkflows/local/input_check'

include { DEPTH_ANALYSIS            as DEPTHANALYSIS        } from '../subworkflows/local/depthanalysis/main'
include { CREATE_PANELS             as CREATEPANELS         } from '../subworkflows/local/createpanels/main'

include { PLOT_DEPTHS               as PLOTDEPTHSALLCONS    } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONS      } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONSCONS  } from '../subworkflows/local/plotdepths/main'

include { MUTATION_PREPROCESSING    as MUT_PREPROCESSING    } from '../subworkflows/local/mutationpreprocessing/main'

include { ENRICHPANELS              as ENRICHPANELS         } from '../subworkflows/local/enrichpanels/main'


include { MUTATION_DENSITY          as MUTDENSITYALL           } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYPROT          } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYNONPROT       } from '../subworkflows/local/mutationdensity/main'
include { MUTATION_DENSITY          as MUTDENSITYSYNONYMOUS    } from '../subworkflows/local/mutationdensity/main'

include { MUTATION_DENSITY          as MUTDENSITYADJUSTED      } from '../subworkflows/local/adjmutdensity/main'

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
include { PLOTTING_QC               as PLOTTINGQC           } from '../subworkflows/local/plotting_qc/main'

include { REGRESSIONS               as REGRESSIONSMUTDENSITY       } from '../subworkflows/local/regressions/main'
include { REGRESSIONS               as REGRESSIONSOMEGA            } from '../subworkflows/local/regressions/main'
include { REGRESSIONS               as REGRESSIONSOMEGAGLOB        } from '../subworkflows/local/regressions/main'

include { DNDS                      as DNDS                 } from '../subworkflows/local/dnds/main'

include { TABIX_BGZIPTABIX_QUERY    as DEPTHSALLCONS        } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSEXONSCONS      } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSPROTCONS       } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSNONPROTCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSINTRONSCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSSYNONYMOUSCONS } from '../modules/nf-core/tabix/bgziptabixquery/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULE
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT CUSTOM MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAF_2_VCF                     as INPUTMAF2VCF                 } from '../modules/local/maf2vcf/main'
include { TABLE_2_GROUP                 as TABLE2GROUP                  } from '../modules/local/table2groups/main'
include { ANNOTATE_DEPTHS               as ANNOTATEDEPTHS               } from '../modules/local/annotatedepth/main'
include { DOWNSAMPLE_DEPTHS             as DOWNSAMPLEDEPTHS             } from '../modules/local/downsample/depths/main'
include { DOWNSAMPLE_DEPTHS             as DOWNSAMPLEDEPTHSALLSAMPLES   } from '../modules/local/downsample/depths/main'

include { TABIX_BGZIPTABIX_QUERY        as QUERYMUTATIONSEXONS          } from '../modules/nf-core/tabix/bgziptabixquery/main'

include { ANALYZE_DEPTHS_GROUPS         as ANALYZEDEPTHSGROUPS          } from '../modules/local/analyzedepths/main'

include { VAF_SMOOTHING                 as VAFSMOOTHING                 } from '../modules/local/vaf_smoothing/main'

include { SELECT_MUTDENSITIES           as SYNMUTDENSITY                } from '../modules/local/select_mutdensity/main'
include { SELECT_MUTDENSITIES           as SYNMUTREADSDENSITY           } from '../modules/local/select_mutdensity/main'
include { SELECT_MUTDENSITIES           as UPDSYNMUTDENSITY             } from '../modules/local/select_mutdensity/main'
include { SELECT_MUTDENSITIES           as UPDSYNMUTREADSDENSITY        } from '../modules/local/select_mutdensity/main'
include { DNDS_PROXY                    as DNDSPROXY                    } from '../modules/local/dnds_proxy/main'

include { MAF_2_VCF                         as MAF2VCF                      } from '../modules/local/maf2vcf/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROMATRIXGENERATOR        } from '../modules/local/signatures/sigprofiler/matrixgenerator/main'
include { SIGPROFILERASSIGNMENT_COSMIC_FIT  as SIGPROFILERASSIGNMENTINDELS  } from '../modules/local/signatures/sigprofiler/assignment/cosmic_fit/main'

include { MUTATIONS_2_SIGNATURES            as MUTS2SIGS                    } from '../modules/local/mutations2sbs/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEEPCSA {

    // // Input channel definitions
    features_table  = channel.fromPath( params.features_table ?: params.input, checkIfExists: true)

    wgs_trinucs     = params.wgs_trinuc_counts
                            ? channel.fromPath( params.wgs_trinuc_counts, checkIfExists: true).first()
                            : channel.empty()
    cosmic_ref      = params.cosmic_ref_signatures
                            ? channel.fromPath( params.cosmic_ref_signatures, checkIfExists: true).first()
                            : channel.empty()
    cosmic_indel_ref   = params.indel_ref_signatures
                            ? channel.fromPath( params.indel_ref_signatures, checkIfExists: true).first()
                            : channel.empty()


    datasets3d      = params.datasets3d
                            ? channel.fromPath( params.datasets3d, checkIfExists: true).first()
                            : channel.empty()
    annotations3d   = params.annotations3d
                            ? channel.fromPath( params.annotations3d, checkIfExists: true).first()
                            : channel.empty()
    seqinfo_df      = params.datasets3d
                            ? channel.fromPath( "${params.datasets3d}/seq_for_mut_prob.tsv", checkIfExists: true).first()
                            : channel.empty()
    cadd_scores     = params.cadd_scores
                            ? channel.of([
                                file(params.cadd_scores, checkIfExists : true),
                                file(params.cadd_scores_ind, checkIfExists : true)
                                ]).first()
                            : channel.empty()

    site_comparison_results         = channel.empty()
    all_compiled_omegas             = channel.empty()
    all_compiled_omegasgloballoc    = channel.value(file("${projectDir}/assets/placeholder_no_file.tsv", checkIfExists: true))
    all_mutdensities_file           = channel.empty()
    all_adjusted_mutdensities_file  = channel.value(file("${projectDir}/assets/placeholder_no_file.tsv", checkIfExists: true))
    all_compiled_stabilities        = channel.empty()

    // if the user wants to use custom gene groups, import the gene groups table
    // otherwise I am using the input csv as a dummy value channel
    custom_groups_table = params.custom_groups_file
                                ? channel.fromPath( params.custom_groups_file, checkIfExists: true).first()
                                : channel.fromPath( params.input )

    // if the user wants to use custom BED file for computing the depths, import the BED file
    // otherwise I am using the input csv as a dummy value channel
    custom_bed_file     = params.custom_bedfile
                                ? channel.fromPath( params.custom_bedfile, checkIfExists: true).first()
                                : channel.fromPath( params.input )

    // Initialize booleans based on user params
    def run_mutabilities    = (params.oncodrivefml || params.oncodriveclustl || params.oncodrive3d)
    def run_mutdensity      = (params.mutationdensity || params.omega)
    def run_profile_all     = (params.profileall || run_mutabilities || run_mutdensity || params.omega)

    // Validate input_maf usage: it requires use_custom_depths to be enabled
    if ( params.input_maf && !params.use_custom_depths ) {
        error "ERROR: '--input_maf' requires '--use_custom_depths true' and a matching '--custom_depths_table'. " +
              "BAM-based depth computation is not supported when providing mutations via a MAF file. " +
              "Please set '--use_custom_depths true' and provide '--custom_depths_table <path>'."
    }

    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    if ( params.input_maf && params.use_custom_depths ) {
        channel.fromPath( params.input_maf, checkIfExists: true)
        .map{ it -> [[id:'input'], it] }
        .first()
        .set { maf_input }
        INPUTMAF2VCF( maf_input )
        INPUTMAF2VCF.out.vcf_files
            .flatten()
            .map { vcf -> 
                def meta = [:]
                meta.id = vcf.baseName
                return [ meta, vcf ]
            }
            .set { sample_vcfs }
        sample_inputs_ch = sample_vcfs
    } else {
        INPUT_CHECK( file(params.input), !params.use_custom_depths )
        sample_inputs_ch = INPUT_CHECK.out.sample_inputs
    }

    // Separate samples and VCFs
    sample_inputs_ch
    .map{ it -> [ "id" : it[0].id ]}
    .set{ meta_samples_alone }

    sample_inputs_ch
    .map{ it -> [it[0], it[1]]}
    .set{ meta_vcfs_alone }


    TABLE2GROUP(features_table)
    grouping_definitions = TABLE2GROUP.out.json_samples.concat(TABLE2GROUP.out.json_groups).concat(TABLE2GROUP.out.json_allgroups).collect()

    // Load group keys from JSON file in 'groups' channel
    TABLE2GROUP.out.json_samples.map { json_path ->
        def json = file(json_path).text
        groovy.json.JsonSlurper.newInstance().parseText(json).keySet()
    }.flatten().unique()
    .set { samples_keys_ch } // this is a channel that contains only the group names as elements of the channel

    TABLE2GROUP.out.json_groups.map { json_path ->
        def json = file(json_path).text
        groovy.json.JsonSlurper.newInstance().parseText(json).keySet()
    }.flatten().unique()
    .set { group_keys_ch } // this is a channel that contains only the group names as elements of the channel


    // Depths and panel creation should be a single subworkflow
    // Depth analysis: compute and plots
    DEPTHANALYSIS(sample_inputs_ch, custom_bed_file)

    // Panels annotation
    CREATEPANELS(DEPTHANALYSIS.out.depths, wgs_trinucs)

    // Mutation preprocessing
    MUT_PREPROCESSING(meta_vcfs_alone,
                        CREATEPANELS.out.all_consensus_bed,
                        CREATEPANELS.out.exons_bed,
                        TABLE2GROUP.out.json_allgroups,
                        group_keys_ch,
                        seqinfo_df,
                        CREATEPANELS.out.added_custom_regions
                        )
    somatic_mutations = MUT_PREPROCESSING.out.somatic_mafs


    positive_selection_results = somatic_mutations

    ANNOTATEDEPTHS(DEPTHANALYSIS.out.depths,
                    CREATEPANELS.out.all_panel,
                    TABLE2GROUP.out.json_allgroups,
                    MUT_PREPROCESSING.out.mask_matrix,
                    file(params.input)
                    )
    ANNOTATEDEPTHS.out.annotated_depths.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ annotated_depths_full }

    // if (params.downsample && params.downsample_proportion < 1) {
    if (params.downsample ){
        DOWNSAMPLEDEPTHS(annotated_depths_full)
        annotated_depths = DOWNSAMPLEDEPTHS.out.downsampled_depths

        DOWNSAMPLEDEPTHSALLSAMPLES(ANNOTATEDEPTHS.out.all_samples_depths)
        all_samples_indv_annotated_depths = DOWNSAMPLEDEPTHSALLSAMPLES.out.downsampled_depths
    } else {
        annotated_depths = annotated_depths_full
        all_samples_indv_annotated_depths = ANNOTATEDEPTHS.out.all_samples_depths
    }

    if (params.plot_depths){
        PLOTDEPTHSALLCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
        PLOTDEPTHSEXONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_bed, CREATEPANELS.out.exons_panel)
    }
    PLOTDEPTHSEXONSCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_consensus_bed, CREATEPANELS.out.exons_consensus_panel)

    // ANALYZEDEPTHSGROUPS should run only when user defines a group list
    if (params.features_groups_list) {
        ANALYZEDEPTHSGROUPS(features_table, PLOTDEPTHSEXONSCONS.out.average_depth_gene_sample, PLOTDEPTHSEXONSCONS.out.average_depth_sample)
    }

    // Enrich regions in consensus panels
    ENRICHPANELS(MUT_PREPROCESSING.out.mutations_all_samples,
                 ANNOTATEDEPTHS.out.all_samples_depths,
                 CREATEPANELS.out.all_consensus_panel,
                 CREATEPANELS.out.nonprot_consensus_panel,
                 CREATEPANELS.out.prot_consensus_panel,
                 CREATEPANELS.out.synonymous_consensus_panel,
                 CREATEPANELS.out.exons_consensus_panel,
                 CREATEPANELS.out.domains_panel_bed, // domains_file
                 CREATEPANELS.out.exons_consensus_bed)


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

    // Intersect BED of all sites with somatic mutations to keep only those mutations in the exons consensus panel
    QUERYMUTATIONSEXONS(somatic_mutations, CREATEPANELS.out.exons_consensus_bed)
    mutations_in_exons = QUERYMUTATIONSEXONS.out.subset


    // Mutational profile
    if ( run_profile_all ){
        MUTPROFILEALL(somatic_mutations, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, wgs_trinucs, TABLE2GROUP.out.json_allgroups)
        all_compiled_stabilities = all_compiled_stabilities.concat(MUTPROFILEALL.out.profile_stabilities.map{ it -> it[1] })
        if (run_mutdensity){
            
            // TODO explore if we could use the ALL consensus panel instead of the exons one
            //      I am suggesting this since we are already distinguishing per consequence type within the script
            MUTDENSITYADJUSTED(somatic_mutations, DEPTHSALLCONS.out.subset, CREATEPANELS.out.exons_consensus_bed, CREATEPANELS.out.exons_consensus_panel, MUTPROFILEALL.out.profile, wgs_trinucs)

            // Concatenate all outputs into a single file
            MUTDENSITYADJUSTED.out.mutdensities.map{ it -> it[1]}.flatten()
            .set{ all_adjusted_mutdensities }
            all_adjusted_mutdensities.collectFile(name: "all_adjusted_mutdensities.tsv", storeDir:"${params.outdir}/mutdensity_adjusted", skip: 1, keepHeader: true).first().set{ all_adjusted_mutdensities_file }

            MUTDENSITYADJUSTED.out.mutdensities_flat.map{ it -> it[1]}.flatten()
            .set{ all_adjusted_mutdensities_flat }
            all_adjusted_mutdensities_flat.collectFile(name: "all_adjusted_mutdensities_flat.tsv", storeDir:"${params.outdir}/mutdensity_adjusted", skip: 1, keepHeader: true)

            channel.of([ [ id: "all_samples" ] ])
            .join( MUTDENSITYADJUSTED.out.mutdensities )
            .set{ all_samples_adj_mutdensity }

            UPDSYNMUTDENSITY(all_samples_adj_mutdensity)

            // UPDSYNMUTREADSDENSITY(all_samples_adj_mutdensity)

            DNDSPROXY(all_adjusted_mutdensities_file, UPDSYNMUTDENSITY.out.mutdensity.first())
        }
    }
    if (params.profilenonprot){
        MUTPROFILENONPROT(somatic_mutations, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, wgs_trinucs, TABLE2GROUP.out.json_allgroups)
        all_compiled_stabilities = all_compiled_stabilities.concat(MUTPROFILENONPROT.out.profile_stabilities.map{ it -> it[1] })
    }
    if (params.profileexons){
        MUTPROFILEEXONS(somatic_mutations, DEPTHSEXONSCONS.out.subset, CREATEPANELS.out.exons_consensus_bed, wgs_trinucs, TABLE2GROUP.out.json_allgroups)
        all_compiled_stabilities = all_compiled_stabilities.concat(MUTPROFILEEXONS.out.profile_stabilities.map{ it -> it[1] })
    }
    if (params.profileintrons){
        DEPTHSINTRONSCONS(annotated_depths, CREATEPANELS.out.introns_consensus_bed)
        MUTPROFILEINTRONS(somatic_mutations, DEPTHSINTRONSCONS.out.subset, CREATEPANELS.out.introns_consensus_bed, wgs_trinucs, TABLE2GROUP.out.json_allgroups)
        all_compiled_stabilities = all_compiled_stabilities.concat(MUTPROFILEINTRONS.out.profile_stabilities.map{ it -> it[1] })
    }
    all_compiled_stabilities.flatten().collectFile(name: "all_profile_stabilities.tsv", storeDir:"${params.outdir}/mutational_profile", skip: 1, keepHeader: true)


    if (run_mutdensity){
        // Mutation Density
        MUTDENSITYALL(somatic_mutations, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, ENRICHPANELS.out.all_consensus_expanded_panel, samples_keys_ch, wgs_trinucs, annotated_depths)
        MUTDENSITYPROT(somatic_mutations, DEPTHSPROTCONS.out.subset, CREATEPANELS.out.prot_consensus_bed, ENRICHPANELS.out.prot_consensus_expanded_panel, samples_keys_ch, wgs_trinucs, annotated_depths)
        MUTDENSITYNONPROT(somatic_mutations, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, ENRICHPANELS.out.nonprot_consensus_expanded_panel, samples_keys_ch, wgs_trinucs, annotated_depths)
        MUTDENSITYSYNONYMOUS(somatic_mutations, DEPTHSSYNONYMOUSCONS.out.subset, CREATEPANELS.out.synonymous_consensus_bed, ENRICHPANELS.out.synonymous_consensus_expanded_panel, samples_keys_ch, wgs_trinucs, annotated_depths)

        // Concatenate all outputs into a single file
        channel.empty()
        .concat(MUTDENSITYALL.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYPROT.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYNONPROT.out.mutdensities.map{ it -> it[1]}.flatten())
        .concat(MUTDENSITYSYNONYMOUS.out.mutdensities.map{ it -> it[1]}.flatten())
        .set{ all_mutdensities }
        all_mutdensities.collectFile(name: "all_mutdensities.tsv", storeDir:"${params.outdir}/mutdensity", skip: 1, keepHeader: true).set{ all_mutdensities_file }

        channel.of([ [ id: "all_samples" ] ])
        .join( MUTDENSITYSYNONYMOUS.out.mutdensities )
        .set{ all_samples_syn_mutdensity }

        SYNMUTDENSITY(all_samples_syn_mutdensity)

        SYNMUTREADSDENSITY(all_samples_syn_mutdensity)

        channel.of([ [ id: "all_samples" ] ])
        .join( somatic_mutations )
        .set{ mutations_all_samples }
        VAFSMOOTHING(mutations_all_samples, all_mutdensities_file, PLOTDEPTHSEXONSCONS.out.average_depth_sample)

    }


    if (run_mutabilities) {
        MUTABILITYALL(mutations_in_exons,
                        annotated_depths,
                        MUTPROFILEALL.out.profile,
                        CREATEPANELS.out.exons_consensus_panel
                        )
        if (params.profilenonprot){
            MUTABILITYNONPROT(mutations_in_exons,
                                annotated_depths,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_panel
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
        ONCODRIVEFMLALL(mutations_in_exons, MUTABILITYALL.out.mutability,
                            CREATEPANELS.out.exons_consensus_panel,
                            cadd_scores, "all"
                        )
        positive_selection_results = positive_selection_results.join(ONCODRIVEFMLALL.out.results_snvs, remainder: true)

        if (params.profilenonprot && params.positive_selection_non_protein_affecting){
            ONCODRIVEFMLNONPROT(mutations_in_exons, MUTABILITYNONPROT.out.mutability,
                                    CREATEPANELS.out.exons_consensus_panel,
                                    cadd_scores, "non_prot_aff"
                                )
        }
    }

    if (params.oncodrive3d){
        // Oncodrive3D
        ONCODRIVE3D(mutations_in_exons, MUTABILITYALL.out.mutability,
                    datasets3d, annotations3d, MUT_PREPROCESSING.out.all_raw_vep_annotation)
        positive_selection_results = positive_selection_results.join(ONCODRIVE3D.out.results, remainder: true)
        positive_selection_results = positive_selection_results.join(ONCODRIVE3D.out.results_pos, remainder: true)
    }

    // if (params.expected_mutated_cells & params.dnds){
    if (params.dnds){
        DNDS(mutations_in_exons,
                    DEPTHSEXONSCONS.out.subset,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel,
                    params.fasta
                    )
    }

    if (params.omega){
        omega_regressions_files = channel.empty()
        omega_regressions_files_gloc = channel.empty()

        // Omega
        OMEGA(mutations_in_exons,
                DEPTHSEXONSCONS.out.subset,
                MUTPROFILEALL.out.profile,
                CREATEPANELS.out.exons_consensus_bed,
                ENRICHPANELS.out.exons_consensus_expanded_panel,
                custom_groups_table,
                SYNMUTDENSITY.out.mutdensity.first(),
                CREATEPANELS.out.panel_annotated_rich,
                "",
                grouping_definitions,
                ENRICHPANELS.out.exons_json_subgenic
                )
        positive_selection_results = positive_selection_results.join(OMEGA.out.results, remainder: true)
        all_compiled_omegas = OMEGA.out.all_compiled
        if (params.omega_mutabilities){
            site_comparison_results = OMEGA.out.site_comparison
        }
        if (params.omega_globalloc){
            positive_selection_results = positive_selection_results.join(OMEGA.out.results_global, remainder: true)
            all_compiled_omegasgloballoc = OMEGA.out.all_globalloc_compiled.first()
        }

        if (params.regressions){
            omega_regressions_files = omega_regressions_files.mix(OMEGA.out.results.map{ it -> it[1] })
            omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGA.out.results_global.map{ it -> it[1] })
        }


        if (params.omega_multi){
              // Omega multi
              OMEGAMULTI(mutations_in_exons,
                          DEPTHSEXONSCONS.out.subset,
                          MUTPROFILEALL.out.profile,
                          CREATEPANELS.out.exons_consensus_bed,
                          ENRICHPANELS.out.exons_consensus_expanded_panel,
                          custom_groups_table,
                          SYNMUTREADSDENSITY.out.mutdensity.first(),
                          CREATEPANELS.out.panel_annotated_rich,
                          ".multi",
                          grouping_definitions,
                          ENRICHPANELS.out.exons_json_subgenic
                          )
              positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results, remainder: true)
              if (params.omega_globalloc){
                  positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results_global, remainder: true)
              }
              if (params.regressions){
                  omega_regressions_files = omega_regressions_files.mix(OMEGAMULTI.out.results.map{ it -> it[1] })
                  omega_regressions_files_gloc = omega_regressions_files_gloc.mix(OMEGAMULTI.out.results_global.map{ it -> it[1] })
              }
        }
        if (params.profilenonprot && params.positive_selection_non_protein_affecting){
            OMEGANONPROT(mutations_in_exons,
                            DEPTHSEXONSCONS.out.subset,
                            MUTPROFILENONPROT.out.profile,
                            CREATEPANELS.out.exons_consensus_bed,
                            ENRICHPANELS.out.exons_consensus_expanded_panel,
                            custom_groups_table,
                            SYNMUTDENSITY.out.mutdensity.first(),
                            CREATEPANELS.out.panel_annotated_rich,
                            ".non_protein_affecting",
                            grouping_definitions,
                            ENRICHPANELS.out.exons_json_subgenic
                            )

            if (params.omega_multi){
                OMEGANONPROTMULTI(mutations_in_exons,
                                    DEPTHSEXONSCONS.out.subset,
                                    MUTPROFILENONPROT.out.profile,
                                    CREATEPANELS.out.exons_consensus_bed,
                                    ENRICHPANELS.out.exons_consensus_expanded_panel,
                                    custom_groups_table,
                                    SYNMUTREADSDENSITY.out.mutdensity.first(),
                                    CREATEPANELS.out.panel_annotated_rich,
                                    ".multi.non_protein_affecting",
                                    grouping_definitions,
                                    ENRICHPANELS.out.exons_json_subgenic
                                    )
            }

        }

    }

    if (params.mutated_cells_vaf && params.omega && params.omega_globalloc){
        MUT_PREPROCESSING.out.somatic_mafs
        .join(meta_samples_alone)
        .set{ sample_mutations_only }

        MUTATEDCELLSVAF(sample_mutations_only,
                                CREATEPANELS.out.exons_consensus_bed,
                                OMEGA.out.results_global,
                                features_table
                                )
    }


    if (params.oncodriveclustl){
        // OncodriveClustl
        if (params.profileall){
            ONCODRIVECLUSTL(somatic_mutations, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_panel)
        }
    }


    if (params.signatures){

        channel.of([ [ id: "all_samples" ] ])
        .join( somatic_mutations )
        .set{ maf2vcf_inputs }

        MAF2VCF(maf2vcf_inputs)
        vcf_files = MAF2VCF.out.vcf_files.flatten().collect()

        SIGPROMATRIXGENERATOR(vcf_files)

        SIGPROMATRIXGENERATOR.out.matrix_ID83
        .map{ it -> [ [id:"all_samples"], "indels", it] }
        .set{ indels_matrix }
        
        SIGPROFILERASSIGNMENTINDELS(indels_matrix, cosmic_indel_ref)

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

    PLOTTINGQC(
                somatic_mutations,
                all_mutdensities_file.first(),
                all_adjusted_mutdensities_file,
                all_compiled_omegasgloballoc,
                PLOTDEPTHSEXONSCONS.out.average_depth_gene_sample.first(),
                all_compiled_omegas,
                // site_comparison_results,
                // ANNOTATEDEPTHS.out.all_samples_depths,
                // TABLE2GROUP.out.json_allgroups,
                CREATEPANELS.out.exons_consensus_panel,
                TABLE2GROUP.out.json_allgroups.first(),
                group_keys_ch
                // CREATEPANELS.out.panel_annotated_rich,
                // seqinfo_df,
                // CREATEPANELS.out.domains_in_panel,
                // DNA2PROTEINMAPPING.out.depths_exons_positions
                )

    if (params.omega || params.oncodrive3d || params.oncodrivefml || params.indels || run_mutdensity) {
        if (params.omega){
            positive_selection_results = positive_selection_results.combine(PLOTTINGQC.out.flagged_omegas)
            if (params.omega_globalloc){
                positive_selection_results = positive_selection_results.combine(all_compiled_omegasgloballoc)
            }
        }
        positive_selection_results_ready = positive_selection_results
        .map { element ->
            def meta = element[0]
            def files = element[1..-1].findAll { it -> it != null }
            [meta, files]
        }
        .filter { _meta, files -> files.size() > 0 }

        PLOTTINGSUMMARY(positive_selection_results_ready,
                        somatic_mutations,
                        all_mutdensities_file.first(),
                        all_adjusted_mutdensities_file,
                        site_comparison_results,
                        ANNOTATEDEPTHS.out.all_samples_depths.first(),
                        TABLE2GROUP.out.json_samples.first(),
                        TABLE2GROUP.out.json_allgroups.first(),

                        CREATEPANELS.out.exons_consensus_panel,
                        ENRICHPANELS.out.exons_consensus_expanded_panel,
                        CREATEPANELS.out.panel_annotated_rich,

                        seqinfo_df,
                        CREATEPANELS.out.domains_in_panel,
                        ENRICHPANELS.out.dna2protein_mapping_depth_exons,
                        group_keys_ch
                        )
    }


    // Regressions
    if (params.regressions){

        // Determine config file based on mode
        def config_file = params.bbgr_mode == 'default' ? null : params.bbgr_custom_config

        if (params.mutationdensity && params.bbgr_mutdensity){
            metric = "mutdensity"
            REGRESSIONSMUTDENSITY(
                config_file ?: params.bbgr_mutdensity_config,
                all_mutdensities_file.first(),
                params.bbgr_metadata,
                params.bbgr_mode,
                metric,
                PLOTTINGQC.out.flagged_omegas.first(),
                TABLE2GROUP.out.json_allgroups.first()
            )
        }

        if (params.omega && params.bbgr_omega){
            metric = "omega"
            REGRESSIONSOMEGA(
                config_file ?: params.bbgr_omega_config,
                OMEGA.out.all_compiled.first(),
                params.bbgr_metadata,
                params.bbgr_mode,
                metric,
                PLOTTINGQC.out.flagged_omegas.first(),
                TABLE2GROUP.out.json_allgroups.first()
            )
        }

        if (params.omega_globalloc && params.bbgr_omegagloballoc){
            metric = "omegagloballoc"
            REGRESSIONSOMEGAGLOB(
                config_file ?: params.bbgr_omegagloballoc_config,
                OMEGA.out.all_globalloc_compiled.first(),
                params.bbgr_metadata,
                params.bbgr_mode,
                metric,
                PLOTTINGQC.out.flagged_omegas.first(),
                TABLE2GROUP.out.json_allgroups.first()
            )
        }
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        channel.topic('versions').unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // Info required for completion email and summary
    def multiqc_report = []

    def summary_params = paramsSummaryMap(workflow)
    workflow_summary    = paramsSummaryMultiqc(summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    methods_description    = methodsDescriptionText(ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    ch_multiqc_files = channel.empty()
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
