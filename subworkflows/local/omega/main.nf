include { SUBSET_MAF                as SUBSETOMEGA                      } from '../../../modules/local/subsetmaf/main'
include { SUBSET_MAF                as SUBSETOMEGAMULTI                 } from '../../../modules/local/subsetmaf/main'

include { TABIX_BGZIPTABIX_QUERY    as QUERYPANEL                       } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { OMEGA_PREPROCESS          as PREPROCESSING                    } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { GROUP_GENES               as GROUPGENES                       } from '../../../modules/local/group_genes/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                        } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIES             } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGA                        } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISON                   } from '../../../modules/local/bbgtools/sitecomparison/main'
include { SITE_COMPARISON           as SITECOMPARISONMULTI              } from '../../../modules/local/bbgtools/sitecomparison/main'
include { PLOT_OMEGASYN_QC          as EVALOMEGAGLOCESTIMATION          } from '../../../modules/local/plot/qc/globalloc_synonymous/main'
include { OMEGA_MULTITEST           as OMEGAMULTIPLETEST                } from '../../../modules/local/omega_multipletesting/main'

include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC           } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC               } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIESGLOBALLOC    } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGAGLOBALLOC               } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOC          } from '../../../modules/local/bbgtools/sitecomparison/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOCMULTI     } from '../../../modules/local/bbgtools/sitecomparison/main'
include { OMEGA_MULTITEST           as OMEGAMULTIPLETESTGLOBALLOC       } from '../../../modules/local/omega_multipletesting/main'
include { HOTSPOTS_SELECTION        as HOTSPOTSSELECTION                } from '../../../modules/local/hotspots_selection/main'
include { HOTSPOTS_SELECTION        as HOTSPOTSSELECTIONGLOBALLOC       } from '../../../modules/local/hotspots_selection/main'

workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    expanded_panel
    custom_gene_groups
    mutationdensities
    panel_captured_rich
    suffix
    grouping_defs
    json_subgenic


    main:

    site_comparison_results = channel.empty()
    global_loc_results      = channel.empty()
    all_gloc_results        = channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYPANEL(panel_captured_rich, bedfile)

    SUBSETOMEGA(mutations)
    SUBSETOMEGAMULTI(mutations)

    SUBSETOMEGA.out.mutations
    .join( depth )
    .join( profile )
    .set{ muts_n_depths_n_profile }

    channel.of([ [ id: "all_samples" ] ])
    .join( profile ).first()
    .set{ all_samples_mut_profile }


    // FIXME here I am using bedfile as a dummy value channel
    PREPROCESSING( muts_n_depths_n_profile,
                    expanded_panel,
                    bedfile,
                    all_samples_mut_profile)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join( depth )
    .set{ preprocess_n_depths }

    GROUPGENES(expanded_panel, custom_gene_groups, json_subgenic)

    ESTIMATOR( preprocess_n_depths, expanded_panel,
                    GROUPGENES.out.json_genes.first(), "${projectDir}/assets/omega_consequences_groupings.json")

    if (params.omega_plot){
        mutations
        .join(ESTIMATOR.out.results)
        .set{mutations_n_omega}

        PLOTOMEGA(mutations_n_omega)
    }

    if (params.omega_mutabilities){
        ABSOLUTEMUTABILITIES(preprocess_n_depths,
                                expanded_panel,
                                GROUPGENES.out.json_genes.first())
        SUBSETOMEGA.out.mutations
        .join(ABSOLUTEMUTABILITIES.out.mutabilities)
        .set{mutations_n_mutabilities}

        SITECOMPARISON(mutations_n_mutabilities,
                        QUERYPANEL.out.subset.first())
        site_comparison_results = SITECOMPARISON.out.comparisons

        SUBSETOMEGAMULTI.out.mutations
        .join(ABSOLUTEMUTABILITIES.out.mutabilities)
        .set{mutations_n_mutabilities_globalloc}

        SITECOMPARISONMULTI(mutations_n_mutabilities_globalloc,
                                QUERYPANEL.out.subset.first())
        // site_comparison_results = site_comparison_results.join(SITECOMPARISONMULTI.out.comparisons, remainder: true)

    }


    if (params.omega_globalloc) {

        PREPROCESSINGGLOBALLOC(muts_n_depths_n_profile,
                                expanded_panel,
                                mutationdensities.first(),
                                all_samples_mut_profile)

        PREPROCESSINGGLOBALLOC.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_globalloc_n_depths }

        ESTIMATORGLOBALLOC(preprocess_globalloc_n_depths,
                            expanded_panel,
                            GROUPGENES.out.json_genes.first(),
                            "${projectDir}/assets/omega_consequences_groupings.json")

        global_loc_results = ESTIMATORGLOBALLOC.out.results
        
        global_loc_results.map{ it -> it[1]}.flatten().set{ all_gloc_indv_results }
        // Keep the concatenated file in the work directory; the corrected output is published below.
        all_gloc_indv_results
        .collectFile(name: "all_omegas${suffix}_global_loc.tsv", skip: 1, keepHeader: true)
        .set{ all_gloc_results_raw }
        OMEGAMULTIPLETESTGLOBALLOC(all_gloc_results_raw, grouping_defs)
        all_gloc_results = OMEGAMULTIPLETESTGLOBALLOC.out.corrected

        PREPROCESSING.out.syn_muts_tsv.map{ it -> it[1]}.flatten().collect().set{ all_syn_muts }
        PREPROCESSINGGLOBALLOC.out.syn_muts_tsv.map{ it -> it[1]}.flatten().collect().set{ all_syn_muts_gloc }
        EVALOMEGAGLOCESTIMATION(all_syn_muts, all_syn_muts_gloc, grouping_defs)

        if (params.omega_plot){
            mutations
            .join(ESTIMATORGLOBALLOC.out.results)
            .set{mutations_n_omegagloloc}

            PLOTOMEGAGLOBALLOC(mutations_n_omegagloloc)
        }

        if (params.omega_mutabilities){
            ABSOLUTEMUTABILITIESGLOBALLOC(preprocess_globalloc_n_depths,
                                            expanded_panel,
                                            GROUPGENES.out.json_genes.first())
            SUBSETOMEGA.out.mutations
            .join(ABSOLUTEMUTABILITIESGLOBALLOC.out.mutabilities)
            .set{mutations_n_mutabilities_globalloc}

            SITECOMPARISONGLOBALLOC(mutations_n_mutabilities_globalloc,
                                    QUERYPANEL.out.subset.first())
            // site_comparison_results = site_comparison_results.join(SITECOMPARISONGLOBALLOC.out.comparisons, remainder: true)


            SUBSETOMEGAMULTI.out.mutations
            .join(ABSOLUTEMUTABILITIESGLOBALLOC.out.mutabilities)
            .set{mutations_n_mutabilities_globalloc}

            SITECOMPARISONGLOBALLOCMULTI(mutations_n_mutabilities_globalloc,
                                            QUERYPANEL.out.subset.first())
            // site_comparison_results = site_comparison_results.join(SITECOMPARISONGLOBALLOCMULTI.out.comparisons, remainder: true)
        }

    }

    site_comparison_results.map { it -> 
        def meta = it[0]
        def all_files = it[1..-1].flatten()
        [meta, all_files]
    }.set{ site_comparison_results_flattened }

    if (params.hotspots_annotation && params.hotspots_definition_file) {
        hotspots_file = channel.fromPath(params.hotspots_definition_file, checkIfExists: true).first()
        
        HOTSPOTSSELECTION(
            site_comparison_results,
            QUERYPANEL.out.subset.first(),
            hotspots_file
        )
        // If needed, we can also collect or emit these results
    }


    ESTIMATOR.out.results.map{ it -> it[1]}.flatten().set{ all_indv_results }
    // Keep the concatenated file in the work directory; the corrected output is published below.
    all_indv_results
    .collectFile(name: "all_omegas${suffix}.tsv", skip: 1, keepHeader: true)
    .set{ all_results_raw }
    OMEGAMULTIPLETEST(all_results_raw, grouping_defs)
    all_results = OMEGAMULTIPLETEST.out.corrected


    emit:
    results                 = ESTIMATOR.out.results
    results_global          = global_loc_results
    expanded_panel          = expanded_panel
    site_comparison         = site_comparison_results_flattened

    all_compiled            = all_results
    all_globalloc_compiled  = all_gloc_results
    // plots = ONCODRIVE3D.out.plots

}
