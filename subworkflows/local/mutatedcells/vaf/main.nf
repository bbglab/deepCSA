include { TABIX_BGZIPTABIX_QUERY as SUBSETMUTATIONS           } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF as SUBSETMUTEPIVAFAM                     } from '../../../../modules/local/subsetmaf/main'
include { MUTATED_GENOMES_FROM_VAF as MUTATEDGENOMESFROMVAFAM } from '../../../../modules/local/mutated_genomes_from_vaf/main'
include { MUTATED_CELLS_FROM_VAF as MUTATEDCELLSFROMVAFAM     } from '../../../../modules/local/mutated_cells_from_vaf/main'


workflow MUTATED_CELLS_VAF {
    take:
    mutations
    bedfile
    omegas
    clinical_features

    main:

    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTEPIVAFAM(SUBSETMUTATIONS.out.subset)

    SUBSETMUTEPIVAFAM.out.mutations
        .join(omegas)
        .set { mutations_n_omega }

    MUTATEDGENOMESFROMVAFAM(mutations_n_omega)
    MUTATEDGENOMESFROMVAFAM.out.mutated_gen_sample.map { it -> it[1] }.collect().map { it -> [[id: "all_samples"], it] }.set { mut_genomes_am_samples }
    MUTATEDCELLSFROMVAFAM(mut_genomes_am_samples, clinical_features)
}
