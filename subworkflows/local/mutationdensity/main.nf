include { TABIX_BGZIPTABIX_QUERY    as QUERYMUTATIONS      } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTDENSITY     } from '../../../modules/local/subsetmaf/main'

include { MUTATION_DENSITY          as MUTDENSITY           } from '../../../modules/local/mut_density/simple/main'


workflow MUTATION_DENSITY{
    take:
    mutations
    depth
    bedfile
    panel

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYMUTATIONS(mutations, bedfile)

    SUBSETMUTDENSITY(QUERYMUTATIONS.out.subset)

    SUBSETMUTDENSITY.out.mutations
    .join(depth)
    .set{ mutations_n_depth }

    MUTDENSITY(mutations_n_depth, panel)


    emit:
    mutdensities = MUTDENSITY.out.mutdensities
}
