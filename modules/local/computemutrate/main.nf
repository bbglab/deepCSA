process MUTRATE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(mutations), path(depth)
    tuple val(meta2), path(consensus_panel)

    output:
    tuple val(meta), path("*.mutrates.tsv"), emit: mutrates
    path  "versions.yml"                   , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
