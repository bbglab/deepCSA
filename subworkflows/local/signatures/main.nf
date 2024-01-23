include { SIGPROFILERASSIGNMENT  } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { MATRIX_CONCAT          } from '../../../modules/local/sig_matrix_concat/main'



workflow SIGNATURES {
    take:
    // inputs
    matrix
    reference_signatures

    main:
    // actual code
    ch_versions = Channel.empty()
    matrix.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ matrix_all_samples }

    MATRIX_CONCAT(matrix_all_samples)

    SIGPROFILERASSIGNMENT(MATRIX_CONCAT.out.wgs_tsv, reference_signatures)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    plots    = SIGPROFILERASSIGNMENT.out.plots  // channel: [ val(meta), file(depths) ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
