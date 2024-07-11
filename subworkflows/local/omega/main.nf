include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'


include { SUBSET_MAF                as SUBSET_OMEGA             } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'

include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'

include { SUBSET_MAF                as SUBSET_OMEGA_EXPANDED    } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGEXP         } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOREXP             } from '../../../modules/local/bbgtools/omega/estimator/main'

include { SUBSET_MAF                as SUBSET_OMEGA_OK          } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGOK          } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOROK              } from '../../../modules/local/bbgtools/omega/estimator/main'

include { SUBSET_MAF                as SUBSET_OMEGA_REDUCED     } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGRED         } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORRED             } from '../../../modules/local/bbgtools/omega/estimator/main'



workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    panel
    globalloc
    vaf_distorsioned

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


    // FIXME: here I am using bedfile as a dummy value channel
    PREPROCESSING( muts_n_depths_n_profile, panel, bedfile)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join( SUBSETDEPTHS.out.subset )
    .set{ preprocess_n_depths }

    ESTIMATOR( preprocess_n_depths, panel)
    ch_versions = ch_versions.mix(ESTIMATOR.out.versions)


    if (globalloc) {
        Channel.of([ [ id: "all_samples" ] ])
        .join( PREPROCESSING.out.syn_muts_tsv )
        .set{ all_samples_syn_muts }


        PREPROCESSINGGLOBALLOC( muts_n_depths_n_profile, panel, all_samples_syn_muts.first() )
        ch_versions = ch_versions.mix(PREPROCESSINGGLOBALLOC.out.versions)

        PREPROCESSINGGLOBALLOC.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_globalloc_n_depths }

        ESTIMATORGLOBALLOC( preprocess_globalloc_n_depths, panel)
        ch_versions = ch_versions.mix(ESTIMATORGLOBALLOC.out.versions)

        global_loc_results = ESTIMATORGLOBALLOC.out.results
    } else {
        global_loc_results = Channel.empty()
    }



    if (vaf_distorsioned) {

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
        PREPROCESSINGEXP( muts_n_depths_n_profile_exp, panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGEXP.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_n_depths_exp }

        ESTIMATOREXP( preprocess_n_depths_exp, panel)
        ch_versions = ch_versions.mix(ESTIMATOREXP.out.versions)




        SUBSET_OMEGA_OK(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSET_OMEGA_OK.out.versions)

        SUBSET_OMEGA_OK.out.mutations
        .join( SUBSETDEPTHS.out.subset )
        .join( profile )
        .set{ muts_n_depths_n_profile_ok }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGOK( muts_n_depths_n_profile_ok, panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGEXP.out.versions)

        PREPROCESSINGOK.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_n_depths_ok }

        ESTIMATOROK( preprocess_n_depths_ok, panel)
        ch_versions = ch_versions.mix(ESTIMATOROK.out.versions)



        SUBSET_OMEGA_REDUCED(all_samples_muts)
        ch_versions = ch_versions.mix(SUBSET_OMEGA_REDUCED.out.versions)

        SUBSET_OMEGA_REDUCED.out.mutations
        .join( SUBSETDEPTHS.out.subset )
        .join( profile )
        .set{ muts_n_depths_n_profile_red }

        // FIXME: here I am using bedfile as a dummy value channel
        PREPROCESSINGRED( muts_n_depths_n_profile_red, panel, bedfile)
        ch_versions = ch_versions.mix(PREPROCESSINGRED.out.versions)

        PREPROCESSINGRED.out.mutabs_n_mutations_tsv
        .join( SUBSETDEPTHS.out.subset )
        .set{ preprocess_n_depths_red }

        ESTIMATORRED( preprocess_n_depths_red, panel)
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
