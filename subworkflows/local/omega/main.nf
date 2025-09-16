include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                as SUBSETOMEGA              } from '../../../modules/local/subsetmaf/main'
include { SUBSET_MAF                as SUBSETOMEGAMULTI         } from '../../../modules/local/subsetmaf/main'



include { TABIX_BGZIPTABIX_QUERY    as SUBSETPANEL              } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { EXPAND_REGIONS            as EXPANDREGIONS            } from '../../../modules/local/expand_regions/main'

include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { GROUP_GENES               as GROUPGENES               } from '../../../modules/local/group_genes/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIES     } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGA                } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISON           } from '../../../modules/local/bbgtools/sitecomparison/main'
include { SITE_COMPARISON           as SITECOMPARISONMULTI      } from '../../../modules/local/bbgtools/sitecomparison/main'
include { PLOT_OMEGASYN_QC          as EVALOMEGAGLOBALLOC       } from '../../../modules/local/plot/qc/globalloc_synonymous/main'

include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIESGLOBALLOC       } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGAGLOBALLOC       } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOC  } from '../../../modules/local/bbgtools/sitecomparison/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOCMULTI  } from '../../../modules/local/bbgtools/sitecomparison/main'
include { PLOT_OMEGASYN_QC          as EVALOMEGAGLOBALLOCMULTI  } from '../../../modules/local/plot/qc/globalloc_synonymous/main'


workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    panel
    custom_gene_groups
    domains_file
    mutationdensities
    complete_panel
    exons_file
    suffix
    grouping_defs


    main:

    // Create a channel for the domains file if omega_autodomains is true
    domains_ch = params.omega_autodomains ? domains_file : []  // .map{ it -> it[1]} : []
    exons_ch = params.omega_autoexons ? exons_file.map{ it -> it[1]} : []

    // Create a channel for the hotspots bedfile if provided
    subgenic_ch = params.omega_subgenic_bedfile ? file(params.omega_subgenic_bedfile) : []


    site_comparison_results = Channel.empty()
    global_loc_results      = Channel.empty()
    all_gloc_results        = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETPANEL(complete_panel, bedfile)

    SUBSETOMEGA(SUBSETMUTATIONS.out.subset)
    SUBSETOMEGAMULTI(SUBSETMUTATIONS.out.subset)

    SUBSETOMEGA.out.mutations
    .join( depth )
    .join( profile )
    .set{ muts_n_depths_n_profile }

    Channel.of([ [ id: "all_samples" ] ])
    .join( profile ).first()
    .set{ all_samples_mut_profile }


    if (params.omega_withingene){
        EXPANDREGIONS(panel, domains_ch, exons_ch, subgenic_ch)
        expanded_panel = EXPANDREGIONS.out.panel_increased.first()
        json_hotspots = EXPANDREGIONS.out.new_regions_json.first()
    } else {
        expanded_panel = panel.first()
        json_hotspots = bedfile.first()
    }

    // FIXME here I am using bedfile as a dummy value channel
    PREPROCESSING( muts_n_depths_n_profile,
                    expanded_panel,
                    bedfile,
                    all_samples_mut_profile)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join( depth )
    .set{ preprocess_n_depths }

    Channel.of([ [ id: "all_samples" ] ])
    .join( PREPROCESSING.out.syn_muts_tsv )
    .set{ all_samples_muts }

    GROUPGENES(all_samples_muts, custom_gene_groups, json_hotspots)

    ESTIMATOR( preprocess_n_depths, expanded_panel, GROUPGENES.out.json_genes.first())

    if (params.omega_plot){
        SUBSETMUTATIONS.out.subset
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
                        SUBSETPANEL.out.subset.first())
        site_comparison_results = SITECOMPARISON.out.comparisons

        SUBSETOMEGAMULTI.out.mutations
        .join(ABSOLUTEMUTABILITIES.out.mutabilities)
        .set{mutations_n_mutabilities_globalloc}

        SITECOMPARISONMULTI(mutations_n_mutabilities_globalloc,
                                SUBSETPANEL.out.subset.first())
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
                            GROUPGENES.out.json_genes.first())

        global_loc_results = ESTIMATORGLOBALLOC.out.results
        
        global_loc_results.map{ it -> it[1]}.flatten().set{ all_gloc_indv_results }
        all_gloc_indv_results.collectFile(name: "all_omegas${suffix}_global_loc.tsv", storeDir:"${params.outdir}/omegagloballoc", skip: 1, keepHeader: true).set{ all_gloc_results }

        PREPROCESSING.out.syn_muts_tsv.map{ it -> it[1]}.flatten().set{ all_syn_muts }
        PREPROCESSINGGLOBALLOC.out.syn_muts_tsv.map{ it -> it[1]}.flatten().set{ all_syn_muts_gloc }
        EVALOMEGAGLOBALLOC(all_syn_muts, all_syn_muts_gloc, grouping_defs)

        if (params.omega_plot){
            SUBSETMUTATIONS.out.subset
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
                                    SUBSETPANEL.out.subset.first())
            // site_comparison_results = site_comparison_results.join(SITECOMPARISONGLOBALLOC.out.comparisons, remainder: true)


            SUBSETOMEGAMULTI.out.mutations
            .join(ABSOLUTEMUTABILITIESGLOBALLOC.out.mutabilities)
            .set{mutations_n_mutabilities_globalloc}

            SITECOMPARISONGLOBALLOCMULTI(mutations_n_mutabilities_globalloc,
                                            SUBSETPANEL.out.subset.first())
            // site_comparison_results = site_comparison_results.join(SITECOMPARISONGLOBALLOCMULTI.out.comparisons, remainder: true)
        }

    }

    site_comparison_results.map {
        def meta = it[0]
        def all_files = it[1..-1].flatten()
        [meta, all_files]
    }.set{ site_comparison_results_flattened }


    ESTIMATOR.out.results.map{ it -> it[1]}.flatten().set{ all_indv_results }
    all_indv_results.collectFile(name: "all_omegas${suffix}.tsv", storeDir:"${params.outdir}/omega", skip: 1, keepHeader: true).set{ all_results }


    emit:
    results                 = ESTIMATOR.out.results
    results_global          = global_loc_results
    expanded_panel          = expanded_panel
    site_comparison         = site_comparison_results_flattened

    all_compiled            = all_results
    all_globalloc_compiled  = all_gloc_results
    // plots = ONCODRIVE3D.out.plots

}
