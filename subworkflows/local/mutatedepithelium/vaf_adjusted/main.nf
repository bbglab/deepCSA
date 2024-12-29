
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
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSETMUTEPIVAF(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSETMUTEPIVAF.out.versions)

    MUTATEDGENOMESFROMVAF(SUBSETMUTEPIVAF.out.mutations)
    ch_versions = ch_versions.mix(MUTATEDGENOMESFROMVAF.out.versions)

    // 1) run driver_priority qmaps

    // 1.1) scripts in driver_priority.py
    // 1.2) qmap generator in prepare_qmaps.ipynb

    // 2) prepare summary tables with prepare_summary_tables.ipynb

    // 3) covered urothelium with TERT mutations: tert_mutations.ipynb

    if (params.all_duplex_counts){
        SUBSETMUTEPIVAFAM(SUBSETMUTATIONS.out.subset)
        ch_versions = ch_versions.mix(SUBSETMUTEPIVAFAM.out.versions)

        MUTATEDGENOMESFROMVAFAM(SUBSETMUTEPIVAFAM.out.mutations)
        ch_versions = ch_versions.mix(MUTATEDGENOMESFROMVAFAM.out.versions)
    }



    emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample
    versions        = ch_versions                // channel: [ versions.yml ]

}
