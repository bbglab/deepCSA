include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { QUERY_BIOMART             as QUERYBIOMART             } from '../../../modules/local/dnds/querybiomart/main'
include { BUILD_REFCDS              as BUILDREFCDS              } from '../../../modules/local/dnds/buildref/main'

include { RUN_DNDS                  as DNDSRUN                  } from '../../../modules/local/dnds/run/main'


workflow DNDS {
    take:
    mutations
    depth
    panel_bedfile
    panel
    fasta
    

    main:

    covariates = params.dnds_covariates ? channel.fromPath( params.dnds_covariates, checkIfExists: true).first() : channel.empty()
    // ref_trans = params.dnds_ref_transcripts ? channel.fromPath( params.dnds_ref_transcripts, checkIfExists: true).first() : channel.empty()

    SUBSET_DNDS(mutations)

    PREPROCESSDEPTHS(depth, panel)

    SUBSET_DNDS.out.mutations
    .join(PREPROCESSDEPTHS.out.depths)
    .set{ mutations_n_depth }

    QUERYBIOMART(panel, panel_bedfile)

    BUILDREFCDS(QUERYBIOMART.out.filtered_biomart, fasta)
    // BUILDREFCDS.out.ref_cds

    DNDSRUN(mutations_n_depth, BUILDREFCDS.out.ref_cds, covariates)

    // // uncomment whenever we can use the custom RefCDS file
    // DNDSRUN(mutations_n_depth, ref_trans, covariates)

    // emit:
    // dnds_values = DNDSRUN.out.dnds_values
}
