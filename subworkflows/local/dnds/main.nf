include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { QUERY_BIOMART             as QUERYBIOMART             } from '../../../modules/local/dnds/querybiomart/main'
include { FILTER_BIOMART            as FILTERBIOMART            } from '../../../modules/local/dnds/filterbiomart/main'
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

    QUERYBIOMART(panel)
    FILTERBIOMART(QUERYBIOMART.out.complete_biomart, panel_bedfile)

    BUILDREFCDS(FILTERBIOMART.out.filtered_biomart, fasta)

    DNDSRUN(mutations_n_depth, BUILDREFCDS.out.ref_cds, covariates)

    DNDSRUN.out.results_cv.map{ it -> it[1]}.flatten().set{ cv_results }
    cv_results.collectFile(name: "all_dNdScv.cv.tsv", storeDir:"${params.outdir}/selection/dndscv/cv", skip: 1, keepHeader: true).set{ all_dndscv_results }
    
    DNDSRUN.out.results_global.map{ it -> it[1]}.flatten().set{ global_results }
    global_results.collectFile(name: "all_dNdScv.global.tsv", storeDir:"${params.outdir}/selection/dndscv/persample", skip: 1, keepHeader: true).set{ all_dndscv_global_results }

    DNDSRUN.out.results_local.map{ it -> it[1]}.flatten().set{ local_results }
    local_results.collectFile(name: "all_dNdScv.local.tsv", storeDir:"${params.outdir}/selection/dndscv/local", skip: 1, keepHeader: true).set{ all_dndscv_local_results }

    emit:
    all_dndscv_results
    all_dndscv_global_results
    all_dndscv_local_results

}
