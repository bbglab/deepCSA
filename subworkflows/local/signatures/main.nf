include { MATRIX_CONCAT                         as MATRIXCONCATWGS          } from '../../../modules/local/sig_matrix_concat/main'
include { SIGPROFILERASSIGNMENT_COSMIC_FIT      as SIGPROFILERASSIGNMENT    } from '../../../modules/local/signatures/sigprofiler/assignment/cosmic_fit/main'
include { SIGPROFILERASSIGNMENT_DECOMPOSE_FIT   as HDPREASSIGNMENT          } from '../../../modules/local/signatures/sigprofiler/assignment/decompose_fit/main'
include { SIGNATURES_PROBABILITIES              as SIGPROBS                 } from '../../../modules/local/combine_sbs/main'
include { HDP_EXTRACTION                        as HDPEXTRACTION            } from '../signatures_hdp/main'


workflow SIGNATURES {
    take:
    matrix_wgs
    reference_signatures
    samples

    main:
    // actual code
    matrix_wgs.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices_wgs }

    MATRIXCONCATWGS(all_matrices_wgs, samples)

    MATRIXCONCATWGS.out.wgs_round_tsv.flatten()
    .map{ it -> 
        def parts = it.name.tokenize('.')
        [ [id: parts[0]], parts[1], it ]
    }
    .set{ named_matrices_wgs_sp }

    MATRIXCONCATWGS.out.wgs_tsv_hdp.flatten()
    .map{ it -> 
        def parts = it.name.tokenize('.')
        [ [id: parts[0]], parts[1], it ]
    }
    .set{ named_matrices_wgs_hdp }

    SIGPROFILERASSIGNMENT(named_matrices_wgs_sp, reference_signatures)

    SIGPROFILERASSIGNMENT.out.mutation_probs.map{ it -> it[1] }.collect().map{ it -> [[ id: "all_samples" ], it]}.set{ all_sigs_probs }
    SIGPROBS(all_sigs_probs)
    SIGPROBS.out.signature_probs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ signature_probs_samples }

    HDPEXTRACTION(named_matrices_wgs_hdp, reference_signatures)

    HDPREASSIGNMENT(named_matrices_wgs_sp, HDPEXTRACTION.out.signatures.first(), reference_signatures)

    emit:
    plots               = SIGPROFILERASSIGNMENT.out.plots       // channel: [ val(meta), file(depths) ]
    mutation_probs      = signature_probs_samples
}