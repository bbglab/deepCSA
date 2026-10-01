include { SUBSET_MAF                as SUBSET_DNDS              } from '../../../modules/local/subsetmaf/main'

include { PREPROCESS_DNDS           as PREPROCESSDEPTHS         } from '../../../modules/local/dnds/preprocess/main'
include { ADAPT_PANEL_REFCDS        as BIOMARTPANEL4REFCDS      } from '../../../modules/local/dnds/adaptpanelrefcds/main'
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
    full_biomart_output = params.dnds_biomart_ref ? channel.fromPath( params.dnds_biomart_ref, checkIfExists: true).first() : channel.empty()

    SUBSET_DNDS(mutations)

    PREPROCESSDEPTHS(depth, panel)

    SUBSET_DNDS.out.mutations
    .join(PREPROCESSDEPTHS.out.depths)
    .set{ mutations_n_depth }

    BIOMARTPANEL4REFCDS(full_biomart_output, panel_bedfile)

    BUILDREFCDS(BIOMARTPANEL4REFCDS.out.filtered_biomart, fasta)

    DNDSRUN(mutations_n_depth, BUILDREFCDS.out.ref_cds, covariates)

    DNDSRUN.out.results_cv.map{ it -> it[1]}.flatten().set{ cv_results }
    cv_results.collectFile(name: "all_dNdScv.cv.tsv", storeDir:"${params.outdir}/selection/dndscv/cv", skip: 1, keepHeader: true).set{ all_dndscv_results }
    
    DNDSRUN.out.results_global.map{ it -> it[1]}.flatten().set{ global_results }
    global_results.collectFile(name: "all_dNdScv.global.tsv", storeDir:"${params.outdir}/selection/dndscv/persample", skip: 1, keepHeader: true).set{ all_dndscv_global_results }

    DNDSRUN.out.results_local.map{ it -> it[1]}.flatten().set{ local_results }
    local_results.collectFile(name: "all_dNdScv.local.tsv", storeDir:"${params.outdir}/selection/dndscv/local", skip: 1, keepHeader: true).set{ all_dndscv_local_results }

    emit:
    dnds_cv_per_sample = DNDSRUN.out.results_cv
    all_dndscv_results
    all_dndscv_global_results
    all_dndscv_local_results

}
