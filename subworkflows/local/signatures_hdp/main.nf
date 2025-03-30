include { SIGPROFILERASSIGNMENT                                     } from '../../../modules/local/signatures/sigprofiler/assignment/main'
include { MATRIX_CONCAT                         as MATRIXCONCATWGS  } from '../../../modules/local/sig_matrix_concat/main'
// include { MATRIX_CONCAT                        as MATRIXCONCAT      } from '../../../modules/local/sig_matrix_concat/main'
include { SIGNATURES_PROBABILITIES              as SIGPROBS         } from '../../../modules/local/combine_sbs/main'
include { MSIGHDP                                                   } from '../../../modules/local/signatures/msighdp/main'
include { RUN_HDP_WRAPPER                       as HDPRUN           } from '../../../modules/local/signatures/hdprun/main'


params.setup_file = null
params.ref_file = null
params.output_base = null
params.norm_file = "NA"
params.prior_file = "NA"
params.n_mut_cutoff = 50
params.sig_activity_threshold = 0
params.cohort_threshold = 0
params.burnin = 10000
params.n_posterior = 100
params.space = 200
params.cpiter = 3
params.iterations = 15
params.cosine_sim_threshold = 0.9
params.max_iter_em = 1000
params.em_frac_threshold = 0.1
params.script_dir = './R/'

process PREPARE_INPUT {

    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/spellegrini87/hdp:1.0.0'

    input:
    tuple val(meta) , path(matrix)
    tuple val(meta2), path(treelayers)

    output:
    tuple val(meta), path("*.hdp.rds"), path("*treelayer.rds")  , emit: input_data
    path "versions.yml"                                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # First, create an R script
    cat <<EOF > process_data.R
    data = read.table("${matrix}", header = FALSE)
    rownames(data) <- data[,c(1)]
    colnames(data) <- data[c(1),]
    data <- data[-c(1),]
    data[,-c(1)] <- sapply(data[,-c(1)], as.numeric)
    data <- data[,-c(1)]
    saveRDS(data, file = "${prefix}.hdp.rds")
    EOF

    # Run the R script
    Rscript process_data.R

    cat <<EOF > process_metadata.R
    data = read.table("${matrix}", header = FALSE)
    data = data[,c(1)]
    data <- data[-c(1),]
    colnames(data) <- c("sample")
    data\$group = "L"
    saveRDS(data, file = "${prefix}.hdp.treelayer.rds")
    EOF

    # Run the R script
    Rscript process_metadata.R
    """
}

process RUN_HDP_CHAIN_SAMPLING {
    input:
    tuple val(meta), path(count_matrix), path(tree_layers)

    output:
    tuple val(meta), path("iteration_dir"), emit: iteration_results

    script:
    """
    Rscript ${params.script_dir}/run_HDP_chainSampling.R \\
        $count_matrix \\
        iteration_dir \\
        $tree_layers \\
        ${params.prior_file} \\
        ${params.n_mut_cutoff} \\
        ${params.burnin} \\
        ${params.n_posterior} \\
        ${params.space} \\
        ${params.cpiter}
    """
}

process PROCESS_HDP_RESULTS {
    input:
    tuple val(meta), path(count_matrix), path(tree_layers)
    path iteration_dir

    output:
    tuple val(meta), path("output_dir"), emit: processed_results

    script:
    """
    Rscript ${params.script_dir}/run_HDP_processing.R \\
        $count_matrix \\
        output_dir \\
        $tree_layers \\
        ${params.prior_file} \\
        ${params.n_mut_cutoff} \\
        ${params.sig_activity_threshold} \\
        ${params.cohort_threshold}
    """
}

process NORMALIZE_SIGNATURES {
    input:
    tuple val(meta), path(output_dir)

    output:
    tuple val(meta), path("normalized_output_dir"), emit: normalized_results

    when:
    params.norm_file != "NA"

    script:
    """
    Rscript ${params.script_dir}/run_HDP_sigNormalising.R \\
        $output_dir \\
        ${params.norm_file}
    cp -r $output_dir normalized_output_dir
    """
}

process COMPARE_NORMALIZED_SIGNATURES {
    input:
    tuple val(meta), path(normalized_output_dir)

    output:
    tuple val(meta), path("compared_normalized_output_dir"), emit: compared_normalized_results

    when:
    params.norm_file != "NA"

    script:
    """
    Rscript ${params.script_dir}/run_HDP_comparing.R \\
        $normalized_output_dir \\
        ${params.cosine_sim_threshold} \\
        ${params.max_iter_em} \\
        ${params.em_frac_threshold} \\
        ${params.ref_file} \\
        "NA" \\
        "TRUE"
    cp -r $normalized_output_dir compared_normalized_output_dir
    """
}

process COMPARE_SIGNATURES {
    input:
    tuple val(meta), path(output_dir)

    output:
    tuple val(meta), path("compared_output_dir"), emit: compared_results

    script:
    """
    Rscript ${params.script_dir}/run_HDP_comparing.R \\
        $output_dir \\
        ${params.cosine_sim_threshold} \\
        ${params.max_iter_em} \\
        ${params.em_frac_threshold} \\
        ${params.ref_file} \\
        "NA" \\
        "FALSE"
    cp -r $output_dir compared_output_dir
    """
}

process COMPARE_CANCER_SIGNATURES {
    input:
    tuple val(meta), path(compared_output_dir)

    output:
    tuple val(meta), path("final_output_dir"), emit: final_results

    script:
    """
    cancer_sigs_file=""  # You might want to parameterize this
    Rscript ${params.script_dir}/run_HDP_comparing.R \\
        $compared_output_dir \\
        ${params.cosine_sim_threshold} \\
        ${params.max_iter_em} \\
        ${params.em_frac_threshold} \\
        ${params.ref_file} \\
        \$cancer_sigs_file \\
        "FALSE"
    cp -r $compared_output_dir final_output_dir
    """
}


workflow HDP_EXTRACTION {

    take:
    matrix
    reference_signatures
    samples

    main:

    input_data_ch = PREPARE_INPUT(params.setup_file)

    iteration_results_ch = RUN_HDP_CHAIN_SAMPLING(input_data_ch)

    processed_results_ch = PROCESS_HDP_RESULTS(input_data_ch, iteration_results_ch)

    if (params.norm_file != "NA") {
        normalized_results_ch = NORMALIZE_SIGNATURES(processed_results_ch)
        compared_normalized_results_ch = COMPARE_NORMALIZED_SIGNATURES(normalized_results_ch)
    }

    compared_results_ch = COMPARE_SIGNATURES(processed_results_ch)

    final_results_ch = COMPARE_CANCER_SIGNATURES(compared_results_ch)


    // actual code
    matrix_wgs.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices_wgs }

    MATRIXCONCATWGS(all_matrices_wgs, samples)

    MATRIXCONCATWGS.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices_wgs }
    MATRIXCONCATWGS.out.wgs_tsv_hdp.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices_wgs_hdp }

    SIGPROFILERASSIGNMENT(named_matrices_wgs, reference_signatures)

    SIGPROFILERASSIGNMENT.out.mutation_probs.map{ it -> it[1] }.collect().map{ it -> [[ id: "all_samples" ], it]}.set{ all_sigs_probs }
    SIGPROBS(all_sigs_probs)
    SIGPROBS.out.signature_probs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ signature_probs_samples }

    HDPRUN(named_matrices_wgs_hdp, named_matrices_wgs, reference_signatures)

    // matrix.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ all_matrices }
    // MATRIXCONCAT(all_matrices, samples)
    // MATRIXCONCAT.out.wgs_tsv.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_matrices }
    // MSIGHDP(matrix)

    emit:
    plots               = SIGPROFILERASSIGNMENT.out.plots       // channel: [ val(meta), file(depths) ]
    // plots_extraction    = MSIGHDP.out.plots                     // channel: [ val(meta), file(depths) ]
    mutation_probs      = signature_probs_samples
}
