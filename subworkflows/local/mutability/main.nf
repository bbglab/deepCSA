
// include { INTERSECT_BED     as BED_INTERSECTALL      } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTPA       } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTNONPA    } from '../../modules/local/bedtools/intersect/main'

include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
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

    // // Intersect BED of all sites with BED of sample filtered sites
    // BED_INTERSECTALL(bedfiles.all_sites, bedfiles.sample_discarded)

    // BED_INTERSECTPA(bedfiles.pa_sites, bedfiles.sample_discarded)

    // BED_INTERSECTNONPA(bedfiles.nonpa_sites, bedfiles.sample_discarded)

    // BED_INTERSECTALL.out.bed
    // .join( BED_INTERSECTPA.out.bed )
    // .join( BED_INTERSECTNONPA.out.bed )
    // .set { all_beds }

    // bedfiles.all_sites
    // .join( bedfiles.pa_sites )
    // .join( bedfiles.nonpa_sites )
    // .set { all_beds }
    
    SUBSETDEPTHS(depth, bedfiles)
    SUBSETMUTATIONS(mutations, bedfiles)

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
