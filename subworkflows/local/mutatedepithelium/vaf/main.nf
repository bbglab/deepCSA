
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS          } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSETMUTEPIVAF          } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAF    } from '../../../../modules/local/mutatedgenomesfromvaf/main'


include { SUBSET_MAF                    as SUBSETMUTEPIVAFAM        } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAFAM  } from '../../../../modules/local/mutatedgenomesfromvaf/main'


workflow MUTATED_EPITHELIUM_VAF {

    take:
    mutations
    bedfile


    main:

    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTEPIVAF(SUBSETMUTATIONS.out.subset)

    MUTATEDGENOMESFROMVAF(SUBSETMUTEPIVAF.out.mutations)

    if (params.all_duplex_counts){
        SUBSETMUTEPIVAFAM(SUBSETMUTATIONS.out.subset)

        MUTATEDGENOMESFROMVAFAM(SUBSETMUTEPIVAFAM.out.mutations)
    }


    // emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample

}
