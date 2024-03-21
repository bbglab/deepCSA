include { CREATECUSTOMBEDFILE           as READSPOSBED          } from '../../../modules/local/createpanels/custombedfile/main'
include { COMPUTE_FRAGMENT_COORDS       as COMPUTEFRAGMENTSCOORDS   } from '../../../modules/local/fragments_from_bam/main'
include { READS_PER_REGION              as READSPERREGION           } from '../../../modules/local/reads_per_region/main'

include { TABIX_BGZIPTABIX_QUERY_INDEX  as SUBSETPILEUP             } from '../../../modules/nf-core/tabix/bgziptabixqueryindex/main'

// include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
// include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
// include { SUBSET_MAF                as SUBSET_MUTRATE           } from '../../../modules/local/subsetmaf/main'
// include { COMPUTE_TRINUCLEOTIDE     as COMPUTETRINUC            } from '../../../modules/local/compute_trinucleotide/main'
// include { COMPUTE_PROFILE           as COMPUTEPROFILE           } from '../../../modules/local/compute_profile/main'




workflow MUTATED_EPITHELIUM {

    take:
    mutations
    bedfile
    panel
    bamfile
    pileup


    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETPILEUP(pileup, bedfile)
    ch_versions = ch_versions.mix(SUBSETPILEUP.out.versions)

    COMPUTEFRAGMENTSCOORDS(bamfile)
    ch_versions = ch_versions.mix(COMPUTEFRAGMENTSCOORDS.out.versions)

    READSPOSBED(panel)
    ch_versions = ch_versions.mix(READSPOSBED.out.versions)

    SUBSETPILEUP.out.subset
    .join(COMPUTEFRAGMENTSCOORDS.out.fragments)
    .set{ pileup_n_fragments }

    READSPERREGION(pileup_n_fragments, READSPOSBED.out.bed)
    ch_versions = ch_versions.mix(READSPERREGION.out.versions)

    emit:
    read_counts = READSPERREGION.out.read_counts
    versions    = ch_versions                // channel: [ versions.yml ]
}
