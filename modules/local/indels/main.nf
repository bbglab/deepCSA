process INDELS_COMPARISON {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations)

    output:
    tuple val(meta), path("*.indels.tsv") , emit: indels
    path  "versions.yml"                  , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    indels_comparison.py \\
                --sample ${prefix} \\
                --filename ${mutations} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta.id}"
    """
    touch ${prefix}.${panel_version}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
