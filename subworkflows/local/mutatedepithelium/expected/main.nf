
include { TABIX_BGZIPTABIX_QUERY        as SUBSETDEPTHS      } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS   } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSET_MUTEPIVAF  } from '../../../../modules/local/subsetmaf/main'
include { EXP_MUTRATE                   as EXPMUTRATE        } from '../../../../modules/local/expected_mutrate/main'




workflow EXPECTED_MUTRATE {

    take:
    mutations
    bedfile
    panel
    depth

    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    // SUBSET_MUTEPIVAF(SUBSETMUTATIONS.out.subset)
    // ch_versions = ch_versions.mix(SUBSET_MUTEPIVAF.out.versions)

    EXPMUTRATE(panel, SUBSETMUTATIONS.out.subset, SUBSETDEPTHS.out.subset)
    ch_versions = ch_versions.mix(EXPMUTRATE.out.versions)

    emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample
    refcds_object   = EXPMUTRATE.out.rds_file.first()
    versions        = ch_versions                // channel: [ versions.yml ]

}
