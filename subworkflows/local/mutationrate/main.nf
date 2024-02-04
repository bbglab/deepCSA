include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_MUTRATE        } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MATRIX            as COMPUTEMATRIX            } from '../../../modules/local/mutation_matrix/main'
include { COMPUTE_TRINUCLEOTIDE     as COMPUTETRINUC            } from '../../../modules/local/compute_trinucleotide/main'

include { COMPUTE_PROFILE           as COMPUTEPROFILE           } from '../../../modules/local/compute_profile/main'


include { MUTRATE as MUTRATE } from '../../../modules/local/computemutrate/main'


workflow MUTATION_RATE{
    take:
    mutations
    depth
    bedfile
    panel

    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSET_MUTRATE(SUBSETMUTATIONS.out.subset)

    SUBSET_MUTRATE.out.mutations
    .join(SUBSETDEPTHS.out.subset)
    .set{ mutations_n_depth }

    MUTRATE(mutations_n_depth, panel)
    // PLOTMUTRATE()

    emit:
    mutrates      = MUTRATE.out.mutrates
    versions = ch_versions                // channel: [ versions.yml ]
}
