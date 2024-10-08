
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_MUTABILITY        } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MUTABILITY        as COMPUTEMUTABILITY        } from '../../../modules/local/compute_mutability/main'

include { TABIX_BGZIPTABIX          as MUTABILITY_BGZIPTABIX    } from '../../../modules/nf-core/tabix/bgziptabix/main'



workflow MUTABILITY {

    take:
    mutations
    depth
    profile
    panel_file
    bedfiles

    main:
    // actual code
    ch_versions = Channel.empty()

    // this line was not needed nor used in the computation of mutabilities,
    // since there is an additional left_join merge in the compute mutabilities code,
    // the depths can be the full ones
    // SUBSETDEPTHS(depth, bedfiles)
    SUBSETMUTATIONS(mutations, bedfiles)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)


    SUBSET_MUTABILITY(SUBSETMUTATIONS.out.subset)
    SUBSET_MUTABILITY.out.mutations
    .join(profile)
    .join(depth)
    .set{ mutations_n_profile_n_depth }

    COMPUTEMUTABILITY( mutations_n_profile_n_depth, panel_file)

    MUTABILITY_BGZIPTABIX( COMPUTEMUTABILITY.out.mutability )


    emit:
    mutability      = MUTABILITY_BGZIPTABIX.out.gz_tbi      // channel: [ val(meta), file(mutabilities), file(mutabilities_index) ]
    versions        = ch_versions                           // channel: [ versions.yml ]
}
