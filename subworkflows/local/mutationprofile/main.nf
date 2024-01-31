
// include { INTERSECT_BED     as BED_INTERSECTALL      } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTPA       } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTNONPA    } from '../../modules/local/bedtools/intersect/main'

include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_MUTPROFILE        } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MATRIX            as COMPUTEMATRIX            } from '../../../modules/local/mutation_matrix/main'
include { COMPUTE_TRINUCLEOTIDE     as COMPUTETRINUC            } from '../../../modules/local/compute_trinucleotide/main'

include { COMPUTE_PROFILE           as COMPUTEPROFILE           } from '../../../modules/local/compute_profile/main'



workflow MUTATIONAL_PROFILE {

    take:
    mutations
    depth
    bedfile

    main:
    // actual code
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSET_MUTPROFILE(SUBSETMUTATIONS.out.subset)

    COMPUTEMATRIX(SUBSET_MUTPROFILE.out.mutations)

    COMPUTEMATRIX.out.per_sample_sigprof
    .join( Channel.of([ [ id: "all_samples" ], [] ]) )
    .map{ it -> [ it[0], it[1]]  }
    .set{ sigprofiler_matrix }

    COMPUTETRINUC(SUBSETDEPTHS.out.subset)

    COMPUTEMATRIX.out.matrix
    .join(COMPUTETRINUC.out.trinucleotides)
    .set{ matrix_n_trinucleotide }

    COMPUTEPROFILE(matrix_n_trinucleotide)

    sigprofiler_empty = Channel.of([])
    sigprofiler_empty
    .concat(COMPUTEPROFILE.out.wgs_sigprofiler)
    .set{ sigprofiler_wgs }


    emit:
    profile         = COMPUTEPROFILE.out.profile            // channel: [ val(meta), file(profile) ]
    matrix_sigprof  = sigprofiler_matrix
    trinucleotides  = COMPUTETRINUC.out.trinucleotides
    wgs_sigprofiler = sigprofiler_wgs
    versions        = ch_versions                           // channel: [ versions.yml ]
}
