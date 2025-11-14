include { PREPARE_INPUT                                         } from '../../../modules/local/signatures/hdp/prepareinput/main'
include { RUN_HDP_CHAIN_SAMPLING                                } from '../../../modules/local/signatures/hdp/chainsampling/main'
include { PROCESS_HDP_RESULTS                                   } from '../../../modules/local/signatures/hdp/process_results/main'
include { COMPARE_SIGNATURES as COMPARESIGNATURES               } from '../../../modules/local/signatures/hdp/compare_sigs/main'




workflow HDP_EXTRACTION {

    take:
    matrix
    reference_signatures

    main:

    Channel.of([ [ id: "samples_matrix" ] ])
    .join( matrix )
    .set{ samples_matrix }

    PREPARE_INPUT(samples_matrix)

    iter_ch = Channel.of(1..15)
    combined_input_ch = PREPARE_INPUT.out.input_data.combine(iter_ch)

    RUN_HDP_CHAIN_SAMPLING(combined_input_ch)

    // Collect all iteration results
    all_iterations_ch = RUN_HDP_CHAIN_SAMPLING.out.iteration_results.groupTuple()

    PROCESS_HDP_RESULTS(PREPARE_INPUT.out.input_data, all_iterations_ch)

    COMPARESIGNATURES(PROCESS_HDP_RESULTS.out.processed_results, reference_signatures)


    // emit:
    // plots               = SIGPROFILERASSIGNMENT.out.plots       // channel: [ val(meta), file(depths) ]
    // // plots_extraction    = MSIGHDP.out.plots                     // channel: [ val(meta), file(depths) ]
    // mutation_probs      = signature_probs_samples
}
