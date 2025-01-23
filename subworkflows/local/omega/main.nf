include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                as SUBSETOMEGA              } from '../../../modules/local/subsetmaf/main'


include { TABIX_BGZIPTABIX_QUERY    as SUBSETPANEL              } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { EXPAND_REGIONS            as EXPANDREGIONS            } from '../../../modules/local/expand_regions/main'

include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { GROUP_GENES               as GROUPGENES               } from '../../../modules/local/group_genes/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIES     } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGA                } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISON           } from '../../../modules/local/bbgtools/sitecomparison/main'

include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_MUTABILITIES        as ABSOLUTEMUTABILITIESGLOBALLOC       } from '../../../modules/local/bbgtools/omega/mutabilities/main'
include { PLOT_OMEGA                as PLOTOMEGAGLOBALLOC       } from '../../../modules/local/plot/omega/main'
include { SITE_COMPARISON           as SITECOMPARISONGLOBALLOC  } from '../../../modules/local/bbgtools/sitecomparison/main'

include { SUBSET_MAF                as SUBSETOMEGA_EXPANDED     } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGEXP         } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOREXP             } from '../../../modules/local/bbgtools/omega/estimator/main'

include { SUBSET_MAF                as SUBSETOMEGA_OK           } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGOK          } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOROK              } from '../../../modules/local/bbgtools/omega/estimator/main'

include { SUBSET_MAF                as SUBSETOMEGA_REDUCED      } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGRED         } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORRED             } from '../../../modules/local/bbgtools/omega/estimator/main'



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
        }

    } else {
        global_loc_results = Channel.empty()
    }



    if (params.omega_vaf_distorsioned) {

        Channel.of([ [ id: "all_samples" ] ])
        .join( SUBSETMUTATIONS.out.subset )
        .set{ all_samples_muts }

        SUBSETOMEGA_EXPANDED(all_samples_muts)

        SUBSETOMEGA_EXPANDED.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_exp }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGEXP( muts_n_depths_n_profile_exp, expanded_panel, bedfile)

        PREPROCESSINGEXP.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_exp }

        ESTIMATOREXP( preprocess_n_depths_exp, expanded_panel, GROUPGENES.out.json_genes.first())




        SUBSETOMEGA_OK(all_samples_muts)

        SUBSETOMEGA_OK.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_ok }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGOK( muts_n_depths_n_profile_ok, expanded_panel, bedfile)

        PREPROCESSINGOK.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_ok }

        ESTIMATOROK( preprocess_n_depths_ok, expanded_panel, GROUPGENES.out.json_genes.first())



        SUBSETOMEGA_REDUCED(all_samples_muts)

        SUBSETOMEGA_REDUCED.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_red }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGRED( muts_n_depths_n_profile_red, expanded_panel, bedfile)

        PREPROCESSINGRED.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_red }

        ESTIMATORRED( preprocess_n_depths_red, expanded_panel, GROUPGENES.out.json_genes.first())

        expanded_results = ESTIMATOREXP.out.results
        concording_results = ESTIMATOROK.out.results
        reduced_results = ESTIMATORRED.out.results
    } else {
        expanded_results = Channel.empty()
    }



    emit:
    results         = ESTIMATOR.out.results
    results_global  = global_loc_results


    // plots = ONCODRIVE3D.out.plots

}
