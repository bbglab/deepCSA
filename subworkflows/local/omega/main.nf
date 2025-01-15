include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'


include { SUBSET_MAF                as SUBSETOMEGA              } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { GROUP_GENES               as GROUPGENES               } from '../../../modules/local/group_genes/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'

include { EXPAND_REGIONS            as EXPANDREGIONS            } from '../../../modules/local/expand_regions/main'
include { PLOT_OMEGA                as PLOTOMEGA                } from '../../../modules/local/plot/omega/main'

include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'
include { PLOT_OMEGA                as PLOTOMEGAGLOBALLOC       } from '../../../modules/local/plot/omega/main'

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
    relative_mutationrates


    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSETOMEGA(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSETOMEGA.out.versions)

    SUBSETOMEGA.out.mutations
    .join( depth )
    .join( profile )
    .set{ muts_n_depths_n_profile }

    Channel.of([ [ id: "all_samples" ] ])
    .join( profile ).first()
    .set{ all_samples_mut_profile }


    if (params.omega_hotspots){
        EXPANDREGIONS(panel, hotspots_file)
        ch_versions = ch_versions.mix(EXPANDREGIONS.out.versions)
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
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join( depth )
    .set{ preprocess_n_depths }

    Channel.of([ [ id: "all_samples" ] ])
    .join( PREPROCESSING.out.syn_muts_tsv )
    .set{ all_samples_muts }

    GROUPGENES(all_samples_muts, custom_gene_groups, json_hotspots)
    ch_versions = ch_versions.mix(GROUPGENES.out.versions)

    ESTIMATOR( preprocess_n_depths, expanded_panel, GROUPGENES.out.json_genes.first())
    ch_versions = ch_versions.mix(ESTIMATOR.out.versions)

    if (params.omega_plot){
        SUBSETMUTATIONS.out.subset
        .join(ESTIMATOR.out.results)
        .set{mutations_n_omega}

        PLOTOMEGA(mutations_n_omega)
        ch_versions = ch_versions.mix(PLOTOMEGA.out.versions)
    }


    if (params.omega_globalloc) {

        PREPROCESSINGGLOBALLOC(muts_n_depths_n_profile,
                                expanded_panel,
                                relative_mutationrates.first(),
                                all_samples_mut_profile)
        ch_versions = ch_versions.mix(PREPROCESSINGGLOBALLOC.out.versions)

        PREPROCESSINGGLOBALLOC.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_globalloc_n_depths }

        ESTIMATORGLOBALLOC(preprocess_globalloc_n_depths,
                            expanded_panel,
                            GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATORGLOBALLOC.out.versions)

        global_loc_results = ESTIMATORGLOBALLOC.out.results

        if (params.omega_plot){
            SUBSETMUTATIONS.out.subset
            .join(ESTIMATORGLOBALLOC.out.results)
            .set{mutations_n_omegagloloc}

            PLOTOMEGAGLOBALLOC(mutations_n_omegagloloc)
            ch_versions = ch_versions.mix(PLOTOMEGAGLOBALLOC.out.versions)
        }

    } else {
        global_loc_results = Channel.empty()
    }



    if (params.omega_vaf_distorsioned) {

        Channel.of([ [ id: "all_samples" ] ])
        .join( SUBSETMUTATIONS.out.subset )
        .set{ all_samples_muts }

        SUBSETOMEGA_EXPANDED(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSETOMEGA_EXPANDED.out.versions)

        SUBSETOMEGA_EXPANDED.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_exp }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGEXP( muts_n_depths_n_profile_exp, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGEXP.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_exp }

        ESTIMATOREXP( preprocess_n_depths_exp, expanded_panel, GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATOREXP.out.versions)




        SUBSETOMEGA_OK(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSETOMEGA_OK.out.versions)

        SUBSETOMEGA_OK.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_ok }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGOK( muts_n_depths_n_profile_ok, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGOK.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_ok }

        ESTIMATOROK( preprocess_n_depths_ok, expanded_panel, GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATOROK.out.versions)



        SUBSETOMEGA_REDUCED(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSETOMEGA_REDUCED.out.versions)

        SUBSETOMEGA_REDUCED.out.mutations
        .join( depth )
        .join( profile )
        .set{ muts_n_depths_n_profile_red }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGRED( muts_n_depths_n_profile_red, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGRED.out.versions)

        PREPROCESSINGRED.out.mutabs_n_mutations_tsv
        .join( depth )
        .set{ preprocess_n_depths_red }

        ESTIMATORRED( preprocess_n_depths_red, expanded_panel, GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATORRED.out.versions)

        expanded_results = ESTIMATOREXP.out.results
        concording_results = ESTIMATOROK.out.results
        reduced_results = ESTIMATORRED.out.results
    } else {
        expanded_results = Channel.empty()
    }



    emit:
    results         = ESTIMATOR.out.results
    results_global  = global_loc_results

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
