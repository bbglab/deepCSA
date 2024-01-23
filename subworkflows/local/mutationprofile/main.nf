
// include { INTERSECT_BED     as BED_INTERSECTALL      } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTPA       } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTNONPA    } from '../../modules/local/bedtools/intersect/main'

include { SUBSET_MAF    as SUBSET_MUTPROFILE } from '../../../modules/local/subsetmaf/main'
include { SUBSET_MAF    as SUBSET_MUTABILITY } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MATRIX          as COMPUTEMATRIX      } from '../../../modules/local/mutation_matrix/main'
include { COMPUTE_PROFILE         as COMPUTEPROFILE     } from '../../../modules/local/compute_profile/main'
include { COMPUTE_TRINUCLEOTIDE   as COMPUTETRINUC      } from '../../../modules/local/compute_trinucleotide/main'
include { COMPUTE_MUTABILITY      as COMPUTEMUTABILITY  } from '../../../modules/local/compute_mutability/main'

include { TABIX_BGZIPTABIX        as MUTABILITY_BGZIPTABIX  } from '../../../modules/nf-core/tabix/bgziptabix/main'



workflow MUTATIONAL_PROFILE {

    take:
    mutations
    depth
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

    SUBSET_MUTPROFILE(mutations)
    SUBSET_MUTABILITY(mutations)

    COMPUTEMATRIX(SUBSET_MUTPROFILE.out.mutations)

    COMPUTEMATRIX.out.per_sample_sigprof
    .join( Channel.of([ [ id: "all_samples" ], [] ]) )
    .map{ it -> [ it[0], it[1]]  }
    .set{ sigprofiler_matrix }


    COMPUTETRINUC(depth)
    COMPUTETRINUC.out.trinucleotides.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_trinucleotides }

    COMPUTEMATRIX.out.matrix
    .join(named_trinucleotides)
    .set{ matrix_n_trinucleotide }

    COMPUTEPROFILE(matrix_n_trinucleotide)

    SUBSET_MUTABILITY.out.mutations
    .join(COMPUTEPROFILE.out.profile)
    .set{ mutations_n_profile }

    COMPUTEMUTABILITY( mutations_n_profile, depth.first() )

    MUTABILITY_BGZIPTABIX( COMPUTEMUTABILITY.out.mutability )

    sigprofiler_empty = Channel.of([])
    sigprofiler_empty
    .concat(COMPUTEPROFILE.out.wgs_sigprofiler)
    .set{ sigprofiler_wgs }


    emit:
    profile         = COMPUTEPROFILE.out.profile            // channel: [ val(meta), file(profile) ]
    mutability      = MUTABILITY_BGZIPTABIX.out.gz_tbi      // channel: [ val(meta), file(mutabilities), file(mutabilities_index) ]
    matrix_sigprof  = sigprofiler_matrix
    trinucleotides  = named_trinucleotides
    wgs_sigprofiler = sigprofiler_wgs
    versions        = ch_versions                           // channel: [ versions.yml ]
}
