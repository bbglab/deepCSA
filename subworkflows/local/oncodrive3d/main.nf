include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF    as SUBSET_ONCODRIVE3D } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVE3D } from '../../../modules/local/bbgtools/oncodrive3d/main'



workflow ONCODRIVE3D_ANALYSIS{

    take:
    mutations
    mutabilities
    bedfile
    datasets

    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_ONCODRIVE3D(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_ONCODRIVE3D.out.versions)

    SUBSET_ONCODRIVE3D.out.mutations
    .join(mutabilities)
    .set{ muts_n_mutability }

    // if (datasets) {
    //      do something
    // } else {
    //     BUILDDATASETS
    // }


    ONCODRIVE3D(muts_n_mutability, datasets)
    ch_versions = ch_versions.mix(ONCODRIVE3D.out.versions)

    emit:
    results     = ONCODRIVE3D.out.csv_genes
    results_pos = ONCODRIVE3D.out.csv_pos

    versions = ch_versions           // channel: [ versions.yml ]

    // plots = ONCODRIVE3D.out.plots

}
