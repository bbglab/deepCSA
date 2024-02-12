include { SUBSET_MAF          as SUBSET_OMEGA  } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS    as PREPROCESSING } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR     as ESTIMATOR     } from '../../../modules/local/bbgtools/omega/estimator/main'



workflow OMEGA_ANALYSIS{

    take:
    muts
    depths
    profile
    annotated_panel

    main:
    ch_versions = Channel.empty()

    SUBSET_OMEGA(muts)
    ch_versions = ch_versions.mix(SUBSET_OMEGA.out.versions)

    SUBSET_OMEGA.out.mutations
    .join( depths )
    .join( profile )
    .set{ muts_n_depths_n_profile }


    PREPROCESSING( muts_n_depths_n_profile, annotated_panel)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    PREPROCESSING.out.mutabs_n_mutations_tsv
    .join(depths)
    .set{ preprocess_n_depths }

    ESTIMATOR( preprocess_n_depths, annotated_panel)
    ch_versions = ch_versions.mix(ESTIMATOR.out.versions)


    emit:
    results     = ESTIMATOR.out.results

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
