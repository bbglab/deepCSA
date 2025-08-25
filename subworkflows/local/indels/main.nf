include { TABIX_BGZIPTABIX_QUERY as SUBSETMUTATIONS } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF as SUBSETINDELS                } from '../../../modules/local/subsetmaf/main'

include { INDELS_COMPARISON as INDELS               } from '../../../modules/local/indels/main'


workflow INDELS_SELECTION {
    take:
    mutations
    bedfile

    main:

    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETINDELS(SUBSETMUTATIONS.out.subset)

    INDELS(SUBSETINDELS.out.mutations)

    emit:
    indels = INDELS.out.indels
}
