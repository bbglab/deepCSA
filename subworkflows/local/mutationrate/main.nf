include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTRATE           } from '../../../modules/local/subsetmaf/main'

include { MUTRATE as MUTRATE } from '../../../modules/local/computemutrate/main'


workflow MUTATION_RATE{
    take:
    mutations
    depth
    bedfile
    panel

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTRATE(SUBSETMUTATIONS.out.subset)

    SUBSETMUTRATE.out.mutations
    .join(depth)
    .set{ mutations_n_depth }

    MUTRATE(mutations_n_depth, panel)

    // PLOTMUTRATE()

    emit:
    mutrates = MUTRATE.out.mutrates
}
