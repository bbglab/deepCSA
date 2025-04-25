process  PROCESS_HDP_RESULTS {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:0.1.0'

    input:
    tuple val(meta), path(count_matrix), path(tree_layers)
    tuple val(meta2), path (iteration_dir)

    output:
    tuple val(meta), path("output_dir"), emit: processed_results
    path "versions.yml"                , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p output_dir/iterations/
    mv hdp_chains_*.RData output_dir/iterations/.
    Rscript /app/HDP_sigExtraction/R/run_HDP_processing.R \\
        $count_matrix \\
        output_dir/ \\
        $tree_layers \\
        ${args}

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

