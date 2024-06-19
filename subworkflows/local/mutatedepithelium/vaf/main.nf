
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS              } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSET_MUTEPIVAF             } from '../../../modules/local/subsetmaf/main'

include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAF        } from '../../../modules/local/mutatedgenomesfromvaf/main'




workflow MUTATED_EPITHELIUM_VAF {

    take:
    mutations
    bedfile


    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_MUTEPIVAF(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_MUTEPIVAF.out.versions)

    MUTATEDGENOMESFROMVAF(SUBSET_MUTEPIVAF.out.mutations)
    ch_versions = ch_versions.mix(MUTATEDGENOMESFROMVAF.out.versions)

    emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample
    versions        = ch_versions                // channel: [ versions.yml ]

}
