include { SUBSET_MAF    as SUBSET_ONCODRIVE3D } from '../../../modules/local/subsetmaf/main'
include { ONCODRIVE3D } from '../../../modules/local/bbgtools/oncodrive3d/main'



workflow ONCODRIVE3D_ANALYSIS{

    take:
    muts

    main:
    ch_versions = Channel.empty()

    //
    // Separate input mutations and mutabilities
    //
    muts.
    map{ it -> [it[0], it[1]]}.
    set{ meta_muts }

    muts.
    map{ it -> [it[0], it[2], it[3]]}.
    set{ meta_mutabs }

    SUBSET_ONCODRIVE3D(meta_muts)
    ch_versions = ch_versions.mix(SUBSET_ONCODRIVE3D.out.versions)

    SUBSET_ONCODRIVE3D.out.mutations
    .join(meta_mutabs)
    .set{ muts_n_mutability }


    ONCODRIVE3D(muts_n_mutability)
    ch_versions = ch_versions.mix(ONCODRIVE3D.out.versions)

    emit:
    results     = ONCODRIVE3D.out.csv_genes
    results_pos = ONCODRIVE3D.out.csv_pos

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
