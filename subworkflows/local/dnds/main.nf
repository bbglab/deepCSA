include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { RUN_DNDS                  as DNDSRUN                  } from '../../../modules/local/dnds/run/main'


workflow DNDS {
    take:
    mutations
    depth
    bedfile
    panel
    covariates
    ref_trans

    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_DNDS(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_DNDS.out.versions)

    PREPROCESSDEPTHS(depth, panel)
    ch_versions = ch_versions.mix(PREPROCESSDEPTHS.out.versions)

    SUBSET_DNDS.out.mutations
    .join(PREPROCESSDEPTHS.out.depths)
    .set{ mutations_n_depth }

    ref_trans.map{it -> [ ["id" : "global_RefCDS"] , it ] }
    .set{ refcds_global }
    DNDSRUN(mutations_n_depth, refcds_global, covariates)
    ch_versions = ch_versions.mix(DNDSRUN.out.versions)

    // // uncomment whenever we can use the custom RefCDS file
    // DNDSRUN(mutations_n_depth, ref_trans, covariates)

    emit:
    // dnds_values = DNDSRUN.out.dnds_values
    versions = ch_versions                // channel: [ versions.yml ]
}
