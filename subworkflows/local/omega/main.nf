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
include { SITE_COMPARISON           as SITECOMPARISONMULTI           } from '../../../modules/local/bbgtools/sitecomparison/main'


include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIESGLOBALLOC       } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGAGLOBALLOC       } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOC  } from '../../../modules/local/bbgtools/sitecomparison/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOCMULTI  } from '../../../modules/local/bbgtools/sitecomparison/main'


workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    panel
    custom_gene_groups
    hotspots_file
    mutationrates
    complete_panel


    main:

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


    if (params.omega_hotspots){
        EXPANDREGIONS(panel, hotspots_file)
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

        SUBSETOMEGAMULTI.out.mutations
        .join(ABSOLUTEMUTABILITIES.out.mutabilities)
        .set{mutations_n_mutabilities_globalloc}

        SITECOMPARISONMULTI(mutations_n_mutabilities_globalloc,
                                SUBSETPANEL.out.subset.first())
    }


    if (params.omega_globalloc) {

        PREPROCESSINGGLOBALLOC(muts_n_depths_n_profile,
                                expanded_panel,
                                mutationrates.first(),
                                all_samples_mut_profile)

        PREPROCESSINGGLOBALLOC.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_globalloc_n_depths }

        ESTIMATORGLOBALLOC(preprocess_globalloc_n_depths,
                            expanded_panel,
                            GROUPGENES.out.json_genes.first())

        global_loc_results = ESTIMATORGLOBALLOC.out.results

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

            SUBSETOMEGAMULTI.out.mutations
            .join(ABSOLUTEMUTABILITIESGLOBALLOC.out.mutabilities)
            .set{mutations_n_mutabilities_globalloc}

            SITECOMPARISONGLOBALLOCMULTI(mutations_n_mutabilities_globalloc,
                                            SUBSETPANEL.out.subset.first())
        }

    } else {
        global_loc_results = Channel.empty()
    }




    emit:
    results         = ESTIMATOR.out.results
    results_global  = global_loc_results


    // plots = ONCODRIVE3D.out.plots

}
