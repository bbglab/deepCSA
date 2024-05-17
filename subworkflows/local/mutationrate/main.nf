include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_MUTRATE           } from '../../../modules/local/subsetmaf/main'

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
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_MUTRATE(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_MUTRATE.out.versions)

    SUBSET_MUTRATE.out.mutations
    .join(SUBSETDEPTHS.out.subset)
    .set{ mutations_n_depth }

    MUTRATE(mutations_n_depth, panel)
    ch_versions = ch_versions.mix(MUTRATE.out.versions)

    // PLOTMUTRATE()

    emit:
    mutrates = MUTRATE.out.mutrates
    versions = ch_versions                // channel: [ versions.yml ]
}
