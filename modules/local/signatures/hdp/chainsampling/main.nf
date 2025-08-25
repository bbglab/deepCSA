process RUN_HDP_CHAIN_SAMPLING {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:0.1.0'

    input:
    tuple val(meta), path(count_matrix), path(tree_layers), val(iter)

    output:
    tuple val(meta), path("iteration_dir/*"), emit: iteration_results
    path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir iteration_dir
    Rscript /app/HDP_sigExtraction/R/run_HDP_chainSampling.R \\
        ${count_matrix} \\
        iteration_dir \\
        ${tree_layers} \\
        ${args} \\
        ${iter}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
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
