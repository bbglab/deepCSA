
process COMPARE_SIGNATURES {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:0.1.0'

    input:
    tuple val(meta), path(output_dir)
    path (reference_signatures)

    output:
    tuple val(meta), path("*.compared_output_dir")  , emit: compared_results
    path "versions.yml"                             , topic: versions


    // when:
    // params.norm_file != "NA"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript /app/HDP_sigExtraction/R/run_HDP_comparing.R \\
        output_dir/ \\
        ${args} \\
        ${reference_signatures} \\
        "NA" \\
        "FALSE"
    cp -r output_dir ${prefix}.compared_output_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

}


// process COMPARE_NORMALIZED_SIGNATURES {
//     input:
//     tuple val(meta), path(normalized_output_dir)

//     output:
//     tuple val(meta), path("compared_normalized_output_dir"), emit: compared_normalized_results

//     when:
//     params.norm_file != "NA"

//     script:
//     """
//     Rscript run_HDP_comparing.R \\
//         $normalized_output_dir \\
//         ${params.cosine_sim_threshold} \\
//         ${params.max_iter_em} \\
//         ${params.em_frac_threshold} \\
//         ${params.ref_file} \\
//         "NA" \\
//         "TRUE"
//     cp -r $normalized_output_dir compared_normalized_output_dir
//     """
// }

// process COMPARE_SIGNATURES {
//     input:
//     tuple val(meta), path(output_dir)

//     output:
//     tuple val(meta), path("compared_output_dir"), emit: compared_results

//     script:
//     """
//     Rscript run_HDP_comparing.R \\
//         $output_dir \\
//         ${params.cosine_sim_threshold} \\
//         ${params.max_iter_em} \\
//         ${params.em_frac_threshold} \\
//         ${params.ref_file} \\
//         "NA" \\
//         "FALSE"
//     cp -r $output_dir compared_output_dir
//     """
// }

// process COMPARE_CANCER_SIGNATURES {
//     input:
//     tuple val(meta), path(compared_output_dir)

//     output:
//     tuple val(meta), path("final_output_dir"), emit: final_results

//     script:
//     """
//     cancer_sigs_file=""  # You might want to parameterize this
//     Rscript run_HDP_comparing.R \\
//         $compared_output_dir \\
//         ${params.cosine_sim_threshold} \\
//         ${params.max_iter_em} \\
//         ${params.em_frac_threshold} \\
//         ${params.ref_file} \\
//         \$cancer_sigs_file \\
//         "FALSE"
//     cp -r $compared_output_dir final_output_dir
//     """
// }
