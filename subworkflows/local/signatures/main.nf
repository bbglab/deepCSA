include { SIGPROFILERASSIGNMENT                                     } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { MATRIX_CONCAT                                             } from '../../../modules/local/sig_matrix_concat/main'
include { SIGNATURES_PROBABILITIES      as SIGPROBS  } from '../../../modules/local/combine_sbs/main'



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

    SIGPROFILERASSIGNMENT.out.mutation_probs.map{ it -> it[1] }.collect().map{ it -> [[ id: "all_samples" ], it]}.set{ all_sigs_probs }
    SIGPROBS(all_sigs_probs)
    SIGPROBS.out.signature_probs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ signature_probs_samples }

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    plots           = SIGPROFILERASSIGNMENT.out.plots               // channel: [ val(meta), file(depths) ]
    mutation_probs  = signature_probs_samples
    versions        = ch_versions                                   // channel: [ versions.yml ]
}
