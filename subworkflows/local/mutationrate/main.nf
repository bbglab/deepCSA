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
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSETMUTRATE(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSETMUTRATE.out.versions)

    SUBSETMUTRATE.out.mutations
    .join(depth)
    .set{ mutations_n_depth }

    MUTRATE(mutations_n_depth, panel)
    ch_versions = ch_versions.mix(MUTRATE.out.versions)

    // PLOTMUTRATE()

    emit:
    mutrates = MUTRATE.out.mutrates
    versions = ch_versions                // channel: [ versions.yml ]
}
