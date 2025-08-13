include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS  } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTDENSITYADJUSTED   } from '../../../modules/local/subsetmaf/main'

include { MUTATION_DENSITY          as MUTDENSITY       } from '../../../modules/local/mut_density/main'


workflow MUTATION_DENSITY{
    take:
    mutations
    depth
    bedfile
    panel
    mutability_file
    trinucleotide_counts

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTDENSITYADJUSTED(SUBSETMUTATIONS.out.subset)

    SUBSETMUTDENSITYADJUSTED.out.mutations
    .join(depth)
    .join(mutability_file)
    .set{ mutations_n_depth_n_profile }

    MUTDENSITY(mutations_n_depth_n_profile, panel, trinucleotide_counts)


    emit:
    mutdensities        = MUTDENSITY.out.mutdensities
    mutdensities_flat   = MUTDENSITY.out.mutdensities_flat
    mutdensities_plots  = MUTDENSITY.out.mutdensities_plots
}
