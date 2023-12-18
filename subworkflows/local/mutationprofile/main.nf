
// include { INTERSECT_BED     as BED_INTERSECTALL      } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTPA       } from '../../modules/local/bedtools/intersect/main'
// include { INTERSECT_BED     as BED_INTERSECTNONPA    } from '../../modules/local/bedtools/intersect/main'

include { COMPUTE_MATRIX     as COMPUTEMATRIX    } from '../../../modules/local/mutation_matrix/main'

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

    COMPUTEMATRIX(mutations)


    emit:
    matrix   = COMPUTEMATRIX.out.matrix   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
