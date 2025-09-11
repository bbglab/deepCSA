
include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTPROFILE         } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MATRIX            as COMPUTEMATRIX            } from '../../../modules/local/mutation_matrix/main'
include { COMPUTE_TRINUCLEOTIDE     as COMPUTETRINUC            } from '../../../modules/local/compute_trinucleotide/main'

include { COMPUTE_PROFILE           as COMPUTEPROFILE           } from '../../../modules/local/compute_profile/main'
include { CONCAT_PROFILES           as CONCATPROFILES           } from '../../../modules/local/concatprofiles/main'


workflow MUTATIONAL_PROFILE {

    take:
    mutations
    depth
    bedfile
    wgs_trinuc

    main:
    // actual code

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTPROFILE(SUBSETMUTATIONS.out.subset)

    COMPUTEMATRIX(SUBSETMUTPROFILE.out.mutations)

    COMPUTEMATRIX.out.per_sample_sigprof
    .join( Channel.of([ [ id: "all_samples" ], [] ]) )
    .map{ it -> [ it[0], it[1]]  }
    .set{ sigprofiler_matrix }

    COMPUTETRINUC(depth)

    COMPUTEMATRIX.out.matrix
    .join(COMPUTETRINUC.out.trinucleotides)
    .set{ matrix_n_trinucleotide }

    COMPUTEPROFILE(matrix_n_trinucleotide, wgs_trinuc)

    sigprofiler_empty = Channel.of([])
    sigprofiler_empty
    .concat(COMPUTEPROFILE.out.wgs_sigprofiler)
    .set{ sigprofiler_wgs }

    compile_all_profiles = COMPUTEPROFILE.out.profile.map{ it -> it[1]}.collect().flatten().map {it -> [ [id:'all_profiles'], it] }
    compile_all_profiles.view()

    CONCATPROFILES(compile_all_profiles)

    emit:
    profile         = COMPUTEPROFILE.out.profile            // channel: [ val(meta), file(profile) ]
    matrix_sigprof  = sigprofiler_matrix
    trinucleotides  = COMPUTETRINUC.out.trinucleotides
    wgs_sigprofiler = sigprofiler_wgs
    compiled_profiles = CONCATPROFILES.out.compiled_profiles
}