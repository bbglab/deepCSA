include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_INDELS           } from '../../../modules/local/subsetmaf/main'

include { INDELS_COMPARISON as INDELS } from '../../../modules/local/indels/main'


workflow INDELS_SELECTION {

    take:
    mutations
    bedfile

    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_INDELS(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_INDELS.out.versions)

    INDELS(SUBSET_INDELS.out.mutations)
    ch_versions = ch_versions.mix(INDELS.out.versions)

    emit:
    indels = INDELS.out.indels
    versions = ch_versions                // channel: [ versions.yml ]
}
