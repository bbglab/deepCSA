include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { RUN_DNDS                  as DNDSRUN                  } from '../../../modules/local/dnds/run/main'


workflow DNDS {
    take:
    mutations
    depth
    bedfile
    panel
    refcds_object
    // covariates

    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_DNDS(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_DNDS.out.versions)

    PREPROCESSDEPTHS(SUBSETDEPTHS.out.subset, panel)
    ch_versions = ch_versions.mix(PREPROCESSDEPTHS.out.versions)

    SUBSET_DNDS.out.mutations
    .join(PREPROCESSDEPTHS.out.depths)
    .set{ mutations_n_depth }

    DNDSRUN(mutations_n_depth, refcds_object)
    ch_versions = ch_versions.mix(DNDSRUN.out.versions)

    // PLOTMUTRATE()

    emit:
    // dnds_values = DNDSRUN.out.dnds_values
    versions = ch_versions                // channel: [ versions.yml ]
}
