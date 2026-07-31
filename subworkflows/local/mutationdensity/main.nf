include { TABIX_BGZIPTABIX_QUERY        as QUERYMUTATIONS       } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                    as SUBSETMUTDENSITY     } from '../../../modules/local/subsetmaf/main'

include { MUTATION_DENSITY              as MUTDENSITY           } from '../../../modules/local/mut_density/simple/main'
include { WG_SCALED_MUTATION_DENSITY    as WGSCALEDMUTDENSITY   } from '../../../modules/local/mut_density/wgscaled/main'


workflow MUTATION_DENSITY{
    take:
    mutations
    depth
    bedfile
    panel
    samples_ch
    wgs_trinucs

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYMUTATIONS(mutations, bedfile)

    SUBSETMUTDENSITY(QUERYMUTATIONS.out.subset)

    SUBSETMUTDENSITY.out.mutations
    .join(depth)
    .set{ mutations_n_depth }

    MUTDENSITY(mutations_n_depth, panel)

    mutations_n_depth
    .map { mut -> tuple(mut[0].id, mut) }
    .join(samples_ch)
    .map { it -> it[1] }
    .set{ mutations_n_depths_samples}
    WGSCALEDMUTDENSITY(mutations_n_depths_samples, bedfile, wgs_trinucs)


    emit:
    mutdensities = MUTDENSITY.out.mutdensities
}
