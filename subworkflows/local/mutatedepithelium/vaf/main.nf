
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS          } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSETMUTEPIVAF          } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAF    } from '../../../../modules/local/mutatedgenomesfromvaf/main'


include { SUBSET_MAF                    as SUBSETMUTEPIVAFAM        } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAFAM  } from '../../../../modules/local/mutatedgenomesfromvaf/main'


workflow MUTATED_EPITHELIUM_VAF {

    take:
    mutations
    bedfile
    omegas    


    main:

    SUBSETMUTATIONS(mutations, bedfile)

    if (params.all_duplex_counts){
        SUBSETMUTEPIVAFAM(SUBSETMUTATIONS.out.subset)

        SUBSETMUTEPIVAFAM.out.mutations
        .join(omegas)
        .set{mutations_n_omega}
        MUTATEDGENOMESFROMVAFAM(mutations_n_omega)
        ch_versions = ch_versions.mix(MUTATEDGENOMESFROMVAFAM.out.versions)
    } else {
        SUBSETMUTEPIVAF(SUBSETMUTATIONS.out.subset)
        ch_versions = ch_versions.mix(SUBSETMUTEPIVAF.out.versions)

        SUBSETMUTEPIVAF.out.mutations
        .join(omegas)
        .set{mutations_n_omega}

        MUTATEDGENOMESFROMVAF(mutations_n_omega)
        ch_versions = ch_versions.mix(MUTATEDGENOMESFROMVAF.out.versions)
    }


    // emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample

}
