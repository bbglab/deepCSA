
include { MATRIX_CONCAT                        as MATRIXCONCATWGS   } from '../../../modules/local/sig_matrix_concat/main'
include { SIGPROFILERASSIGNMENT                                     } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { SIGNATURES_PROBABILITIES              as SIGPROBS         } from '../../../modules/local/combine_sbs/main'
// include { MSIGHDP                                                   } from '../../../modules/local/signatures/msighdp/main'
include { HDP_EXTRACTION                        as HDPEXTRACTION    } from '../signatures_hdp/main'


workflow SIGNATURES {
    take:
    matrix_wgs
    reference_signatures
    samples

    main:
    // actual code
    matrix_wgs.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices_wgs }

    MATRIXCONCATWGS(all_matrices_wgs, samples)

    MATRIXCONCATWGS.out.wgs_tsv.flatten()
    .map{ it -> 
        def parts = it.name.tokenize('.')
        [ [id: parts[0]], parts[1], it ]
    }
    .set{ named_matrices_wgs }

    MATRIXCONCATWGS.out.wgs_tsv_hdp.flatten()
    .map{ it -> 
        def parts = it.name.tokenize('.')
        [ [id: parts[0]], parts[1], it ]
    }
    .set{ named_matrices_wgs_hdp }

    SIGPROFILERASSIGNMENT(named_matrices_wgs, reference_signatures)

    SIGPROFILERASSIGNMENT.out.mutation_probs.map{ it -> it[1] }.collect().map{ it -> [[ id: "all_samples" ], it]}.set{ all_sigs_probs }
    SIGPROBS(all_sigs_probs)
    SIGPROBS.out.signature_probs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ signature_probs_samples }

    HDPEXTRACTION(named_matrices_wgs_hdp, reference_signatures)

    // matrix.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices }
    // MATRIXCONCAT(all_matrices, samples)
    // MATRIXCONCAT.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices }
    // MSIGHDP(matrix)

    emit:
    plots               = SIGPROFILERASSIGNMENT.out.plots       // channel: [ val(meta), file(depths) ]
    // plots_extraction    = MSIGHDP.out.plots                     // channel: [ val(meta), file(depths) ]
    mutation_probs      = signature_probs_samples
}
