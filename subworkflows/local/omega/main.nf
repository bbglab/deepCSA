include { SUBSET_MAF          as SUBSET_OMEGA  } from '../../../modules/local/subsetmaf/main'
include { OMEGA_PREPROCESS    as PREPROCESSING } from '../../../modules/local/bbgtools/omega/preprocess/main'
include { OMEGA_ESTIMATOR     as ESTIMATOR     } from '../../../modules/local/bbgtools/omega/estimator/main'



workflow OMEGA_ANALYSIS{

    take:
    muts
    annotated_panel

    main:
    ch_versions = Channel.empty()

    //
    // Separate input mutations and mutabilities
    //
    muts.
    map{ it -> [it[0], it[1]]}.
    set{ meta_muts }

    muts.
    map{ it -> [it[0], it[2]]}.
    set{ meta_profile }

    muts.
    map{ it -> [it[0], it[2], it[3]]}.
    set{ meta_profile_n_depths }


    SUBSET_OMEGA(meta_muts)
    ch_versions = ch_versions.mix(SUBSET_OMEGA.out.versions)

    SUBSET_OMEGA.out.mutations
    .join(meta_profile_n_depths)
    .set{ muts_n_profile_n_depths }

    // TODO
    // Here we would need to add a module or something that puts together all the mutational profiles of the different samples
    // and the same for the depths and the mutations all into the same dataframe
    // this would be a value channel with 3 elements, mutations, depths and mutation profile


    PREPROCESSING( muts_n_profile_n_depths, annotated_panel)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    ESTIMATOR( PREPROCESSING.out.mutabs_n_mutations_tsv, annotated_panel)
    ch_versions = ch_versions.mix(ESTIMATOR.out.versions)


    emit:
    results     = ESTIMATOR.out.results

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
