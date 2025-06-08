process MUTRATE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations), path(depth)
    tuple val(meta2), path(consensus_panel)

    output:
    tuple val(meta), path("*.mutrates.tsv"), emit: mutrates
    path  "versions.yml"                   , topic: versions



    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    compute_mutrate.py \\
                ${mutations} \\
                ${depth} \\
                ${consensus_panel} \\
                ${prefix} \\
                ${panel_version};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
