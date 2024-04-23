include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'


include { SUBSET_MAF                as SUBSET_OMEGA             } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS          as PREPROCESSING            } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATOR                } from '../../../modules/local/bbgtools/omega/estimator/main'
include { OMEGA_PREPROCESS          as PREPROCESSINGGLOBALLOC   } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR           as ESTIMATORGLOBALLOC       } from '../../../modules/local/bbgtools/omega/estimator/main'


workflow OMEGA_ANALYSIS{

    take:
    mutations
    depth
    profile
    bedfile
    panel
    globalloc

    main:
    ch_versions = Channel.empty()

    // TODO
    // We should revise how is tabix doing the subset
    //      whether it include the ends or not
    //      and see how to take this into account for the other analysis

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



    emit:
    results         = ESTIMATOR.out.results
    results_global  = global_loc_results

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
