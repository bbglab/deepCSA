include { SIGPROFILERASSIGNMENT                                     } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { MATRIX_CONCAT                        as MATRIXCONCATWGS   } from '../../../modules/local/sig_matrix_concat/main'
// include { MATRIX_CONCAT                        as MATRIXCONCAT      } from '../../../modules/local/sig_matrix_concat/main'
include { SIGNATURES_PROBABILITIES             as SIGPROBS     } from '../../../modules/local/combine_sbs/main'
include { MSIGHDP                                                   } from '../../../modules/local/signatures/msighdp/main'



workflow SIGNATURES {
    take:
    matrix
    matrix_wgs
    reference_signatures
    samples

    main:
    // actual code
    ch_versions = Channel.empty()
    matrix_wgs.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices_wgs }

    MATRIXCONCATWGS(all_matrices_wgs, samples)
    ch_versions = ch_versions.mix(MATRIXCONCATWGS.out.versions)

    MATRIXCONCATWGS.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices_wgs }

    SIGPROFILERASSIGNMENT(named_matrices_wgs, reference_signatures)
    ch_versions = ch_versions.mix(SIGPROFILERASSIGNMENT.out.versions)

    SIGPROFILERASSIGNMENT.out.mutation_probs.map{ it -> it[1] }.collect().map{ it -> [[ id: "all_samples" ], it]}.set{ all_sigs_probs }
    SIGPROBS(all_sigs_probs)
    ch_versions = ch_versions.mix(SIGPROBS.out.versions)
    SIGPROBS.out.signature_probs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ signature_probs_samples }


    // matrix.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices }
    // MATRIXCONCAT(all_matrices, samples)
    // MATRIXCONCAT.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices }
    MSIGHDP(matrix)
    ch_versions = ch_versions.mix(MSIGHDP.out.versions)

    emit:
    plots               = SIGPROFILERASSIGNMENT.out.plots               // channel: [ val(meta), file(depths) ]
    plots_extraction    = MSIGHDP.out.plots               // channel: [ val(meta), file(depths) ]
    mutation_probs      = signature_probs_samples
    versions            = ch_versions                                   // channel: [ versions.yml ]
}
