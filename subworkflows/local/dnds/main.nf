include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { QUERY_BIOMART             as QUERYBIOMART             } from '../../../modules/local/dnds/querybiomart/main'
include { BUILD_REFCDS              as BUILDREFCDS              } from '../../../modules/local/dnds/buildref/main'

include { RUN_DNDS                  as DNDSRUN                  } from '../../../modules/local/dnds/run/main'


workflow DNDS {
    take:
    mutations
    depth
    bedfile
    panel
    covariates
    fasta

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


    QUERYBIOMART(panel, bedfile)

    BUILDREFCDS(QUERYBIOMART.out.filtered_biomart, fasta)
    // BUILDREFCDS.out.ref_cds

    DNDSRUN(mutations_n_depth, BUILDREFCDS.out.ref_cds, covariates)
    ch_versions = ch_versions.mix(DNDSRUN.out.versions)

    // // uncomment whenever we can use the custom RefCDS file
    // DNDSRUN(mutations_n_depth, ref_trans, covariates)

    emit:
    // dnds_values = DNDSRUN.out.dnds_values
    versions = ch_versions                // channel: [ versions.yml ]
}
