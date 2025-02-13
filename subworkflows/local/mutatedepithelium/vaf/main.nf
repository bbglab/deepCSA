
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS          } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSETMUTEPIVAF          } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAF    } from '../../../../modules/local/mutatedgenomesfromvaf/main'
include { MUTATED_CELLS_FROM_VAF        as MUTATEDCELLSFROMVAF      } from '../../../../modules/local/mutatedcellsfromvaf/main'

include { SUBSET_MAF                    as SUBSETMUTEPIVAFAM        } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF      as MUTATEDGENOMESFROMVAFAM  } from '../../../../modules/local/mutatedgenomesfromvaf/main'
include { MUTATED_CELLS_FROM_VAF        as MUTATEDCELLSFROMVAFAM    } from '../../../../modules/local/mutatedcellsfromvaf/main'


workflow MUTATED_EPITHELIUM_VAF {

    take:
    mutations
    bedfile
    omegas
    clinical_features


    main:

    SUBSETMUTATIONS(mutations, bedfile)

    if (params.all_duplex_counts){
        SUBSETMUTEPIVAFAM(SUBSETMUTATIONS.out.subset)

        SUBSETMUTEPIVAFAM.out.mutations
        .join(omegas)
        .set{mutations_n_omega}

        MUTATEDGENOMESFROMVAFAM(mutations_n_omega)
        MUTATEDGENOMESFROMVAFAM.out.mutated_gen_sample.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ mut_genomes_am_samples }
        MUTATEDCELLSFROMVAFAM(mut_genomes_am_samples, clinical_features)

    } else {
        SUBSETMUTEPIVAF(SUBSETMUTATIONS.out.subset)

        SUBSETMUTEPIVAF.out.mutations
        .join(omegas)
        .set{mutations_n_omega}

        MUTATEDGENOMESFROMVAF(mutations_n_omega)
        MUTATEDGENOMESFROMVAF.out.mutated_gen_sample.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ mut_genomes_samples }
        MUTATEDCELLSFROMVAF(mut_genomes_samples, clinical_features)

    }




    // emit:
    // TODO add some other output
    // mut_epi_sample  = MUTATEDCELLSFROMVAFAM.out.mutated_cells_sample

}
