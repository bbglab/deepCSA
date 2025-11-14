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

    main:

    covariates = params.dnds_covariates ? Channel.fromPath( params.dnds_covariates, checkIfExists: true).first() : Channel.empty()
    ref_trans = params.dnds_ref_transcripts ? Channel.fromPath( params.dnds_ref_transcripts, checkIfExists: true).first() : Channel.empty()


    SUBSETMUTATIONS(mutations, bedfile)

    SUBSET_DNDS(SUBSETMUTATIONS.out.subset)

    PREPROCESSDEPTHS(depth, panel)

    SUBSET_DNDS.out.mutations
    .join(PREPROCESSDEPTHS.out.depths)
    .set{ mutations_n_depth }

    ref_trans.map{it -> [ ["id" : "global_RefCDS"] , it ] }
    .set{ refcds_global }
    DNDSRUN(mutations_n_depth, refcds_global, covariates)

    // // uncomment whenever we can use the custom RefCDS file
    // DNDSRUN(mutations_n_depth, ref_trans, covariates)

    // emit:
    // dnds_values = DNDSRUN.out.dnds_values
}
