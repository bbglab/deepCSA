include { SIGPROFILERASSIGNMENT  } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { MATRIX_CONCAT          } from '../../../modules/local/sig_matrix_concat/main'



workflow SIGNATURES {
    take:
    matrix
    reference_signatures
    samples

    main:
    // actual code
    ch_versions = Channel.empty()
    matrix.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices }

    MATRIX_CONCAT(all_matrices, samples)

    MATRIX_CONCAT.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices }

    SIGPROFILERASSIGNMENT(named_matrices, reference_signatures)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    plots           = SIGPROFILERASSIGNMENT.out.plots               // channel: [ val(meta), file(depths) ]
    mutation_probs  = SIGPROFILERASSIGNMENT.out.mutation_probs
    versions        = ch_versions                                   // channel: [ versions.yml ]
}
