
include { TABIX_BGZIPTABIX_QUERY    as QUERYMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTPROFILE         } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MATRIX            as COMPUTEMATRIX            } from '../../../modules/local/mutation_matrix/main'
include { COMPUTE_TRINUCLEOTIDE     as COMPUTETRINUC            } from '../../../modules/local/compute_trinucleotide/main'

include { COMPUTE_PROFILE           as COMPUTEPROFILE           } from '../../../modules/local/compute_profile/main'
include { COMPUTE_PROFILE           as COMPUTEPROFILECOHORT     } from '../../../modules/local/compute_profile/main'
include { CONCAT_PROFILES           as CONCATPROFILES           } from '../../../modules/local/concatprofiles/main'


workflow MUTATIONAL_PROFILE {

    take:
    mutations
    depth
    bedfile
    wgs_trinuc
    all_groups

    main:
    // actual code

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYMUTATIONS(mutations, bedfile)

    SUBSETMUTPROFILE(QUERYMUTATIONS.out.subset)

    COMPUTEMATRIX(SUBSETMUTPROFILE.out.mutations)

    COMPUTEMATRIX.out.per_sample_sigprof
    .join( channel.of([ [ id: "all_samples" ], [] ]) )
    .map{ it -> [ it[0], it[1]]  }
    .set{ sigprofiler_matrix }

    COMPUTETRINUC(depth, wgs_trinuc)

    COMPUTEMATRIX.out.matrix
    .join(COMPUTETRINUC.out.trinucleotides)
    .set{ matrix_n_trinucleotide }

    dummy_file = wgs_trinuc.map{ it -> [ [ id: "dummy_file" ], it ]  }

    if (params.profile_smoothing) {
        matrix_n_trinucleotide
        .join( channel.of([ [ id: "all_samples" ] ]) )
        .set{ matrix_n_trinucleotide_all }

        COMPUTEPROFILECOHORT(matrix_n_trinucleotide_all, wgs_trinuc, dummy_file)
        cohort_profile = COMPUTEPROFILECOHORT.out.wgs_proportions.first()
    } else {
        cohort_profile = dummy_file
    }
    COMPUTEPROFILE(matrix_n_trinucleotide, wgs_trinuc, cohort_profile)

    sigprofiler_empty = channel.of([])
    sigprofiler_empty
    .concat(COMPUTEPROFILE.out.wgs_sigprofiler)
    .set{ sigprofiler_wgs }

    compile_all_profiles = COMPUTEPROFILE.out.profile.map{ it -> it[1] }.collect().map { files -> [ [id:'all_samples'], files ] }
    CONCATPROFILES(compile_all_profiles, all_groups)

    compile_stabilities = COMPUTEPROFILE.out.profile_stability.map{ it -> it[1] }.collect().map { files -> [ [id:'all_stabilities'], files ] }


    emit:
    profile             = COMPUTEPROFILE.out.profile            // channel: [ val(meta), file(profile) ]
    matrix_sigprof      = sigprofiler_matrix
    trinucleotides      = COMPUTETRINUC.out.trinucleotides
    wgs_sigprofiler     = sigprofiler_wgs
    compiled_profiles   = CONCATPROFILES.out.compiled_profiles
    profile_stabilities = compile_stabilities
}