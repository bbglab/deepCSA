process RUN_HDP_CHAIN_SAMPLING {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:latest'


    input:
    tuple val(meta), path(count_matrix), path(tree_layers), val(iter)

    output:
    tuple val(meta), path("iteration_dir/*"), emit: iteration_results
    path "versions.yml"                     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript run_HDP_chainSampling.R \\
        $count_matrix \\
        iteration_dir \\
        $tree_layers \\
        ${params.prior_file} \\
        ${params.n_mut_cutoff} \\
        ${params.burnin} \\
        ${params.n_posterior} \\
        ${params.space} \\
        ${params.cpiter} \\
        ${iter}
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

