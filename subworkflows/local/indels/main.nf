include { TABIX_BGZIPTABIX_QUERY    as QUERYMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETINDELS           } from '../../../modules/local/subsetmaf/main'

include { INDELS_COMPARISON as INDELS } from '../../../modules/local/indels/main'


workflow INDELS_SELECTION {

    take:
    mutations
    bedfile

    main:

    QUERYMUTATIONS(mutations, bedfile)

    SUBSETINDELS(QUERYMUTATIONS.out.subset)

    INDELS(SUBSETINDELS.out.mutations)

    emit:
    indels = INDELS.out.indels
}
