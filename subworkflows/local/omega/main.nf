include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'


include { SUBSET_MAF                as SUBSET_OMEGA             } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { GROUP_GENES               as GROUPGENES               } from '../../../modules/local/group_genes/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'

if (params.omega_hotspots){
    include { EXPAND_REGIONS           as EXPANDREGIONS                } from '../../../modules/local/expand_regions/main'
}


if (params.omega_plot){
    include { PLOT_OMEGA           as PLOTOMEGA                     } from '../../../modules/local/plot/omega/main'
}

if (params.omega_globalloc){
    include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
    include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'
    if (params.omega_plot){
        include { PLOT_OMEGA           as PLOTOMEGAGLOBALLOC            } from '../../../modules/local/plot/omega/main'
    }
}

if (params.omega_vaf_distorsioned){
    include { SUBSET_MAF                as SUBSET_OMEGA_EXPANDED    } from '../../../modules/local/subsetmaf/main'
    include { OMEGA_PREPROCESS          as PREPROCESSINGEXP         } from '../../../modules/local/bbgtools/omega/preprocess/main'
    include { OMEGA_ESTIMATOR           as ESTIMATOREXP             } from '../../../modules/local/bbgtools/omega/estimator/main'

    include { SUBSET_MAF                as SUBSET_OMEGA_OK          } from '../../../modules/local/subsetmaf/main'
    include { OMEGA_PREPROCESS          as PREPROCESSINGOK          } from '../../../modules/local/bbgtools/omega/preprocess/main'
    include { OMEGA_ESTIMATOR           as ESTIMATOROK              } from '../../../modules/local/bbgtools/omega/estimator/main'

    include { SUBSET_MAF                as SUBSET_OMEGA_REDUCED     } from '../../../modules/local/subsetmaf/main'
    include { OMEGA_PREPROCESS          as PREPROCESSINGRED         } from '../../../modules/local/bbgtools/omega/preprocess/main'
    include { OMEGA_ESTIMATOR           as ESTIMATORRED             } from '../../../modules/local/bbgtools/omega/estimator/main'
}





workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    panel
    custom_gene_groups
    hotspots_file


    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_OMEGA(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_OMEGA.out.versions)

    SUBSET_OMEGA.out.mutations
    .join( SUBSETDEPTHS.out.subset )
    .join( profile )
    .set{ muts_n_depths_n_profile }


    if (params.omega_hotspots){
        EXPANDREGIONS(panel, hotspots_file)
        ch_versions = ch_versions.mix(EXPANDREGIONS.out.versions)
        expanded_panel = EXPANDREGIONS.out.panel_increased
        json_hotspots = EXPANDREGIONS.out.new_regions_json
    } else {
        expanded_panel = panel
        json_hotspots = bedfile
    }

    // FIXME: here I am using bedfile as a dummy value channel
    PREPROCESSING( muts_n_depths_n_profile, expanded_panel, bedfile)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join( SUBSETDEPTHS.out.subset )
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
        Channel.of([ [ id: "all_samples" ] ])
        .join( PREPROCESSING.out.syn_muts_tsv )
        .set{ all_samples_syn_muts }


        PREPROCESSINGGLOBALLOC( muts_n_depths_n_profile, expanded_panel, all_samples_syn_muts.first() )
        ch_versions = ch_versions.mix(PREPROCESSINGGLOBALLOC.out.versions)

        PREPROCESSINGGLOBALLOC.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_globalloc_n_depths }

        ESTIMATORGLOBALLOC( preprocess_globalloc_n_depths, expanded_panel, GROUPGENES.out.json_genes.first())
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

        SUBSET_OMEGA_EXPANDED(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSET_OMEGA_EXPANDED.out.versions)

        SUBSET_OMEGA_EXPANDED.out.mutations
        .join( SUBSETDEPTHS.out.subset )
        .join( profile )
        .set{ muts_n_depths_n_profile_exp }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGEXP( muts_n_depths_n_profile_exp, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGEXP.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_n_depths_exp }

        ESTIMATOREXP( preprocess_n_depths_exp, expanded_panel, GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATOREXP.out.versions)




        SUBSET_OMEGA_OK(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSET_OMEGA_OK.out.versions)

        SUBSET_OMEGA_OK.out.mutations
        .join( SUBSETDEPTHS.out.subset )
        .join( profile )
        .set{ muts_n_depths_n_profile_ok }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGOK( muts_n_depths_n_profile_ok, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGOK.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_n_depths_ok }

        ESTIMATOROK( preprocess_n_depths_ok, expanded_panel, GROUPGENES.out.json_genes.first())
        ch_versions = ch_versions.mix(ESTIMATOROK.out.versions)



        SUBSET_OMEGA_REDUCED(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSET_OMEGA_REDUCED.out.versions)

        SUBSET_OMEGA_REDUCED.out.mutations
        .join( SUBSETDEPTHS.out.subset )
        .join( profile )
        .set{ muts_n_depths_n_profile_red }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGRED( muts_n_depths_n_profile_red, expanded_panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGRED.out.versions)

        PREPROCESSINGRED.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
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
