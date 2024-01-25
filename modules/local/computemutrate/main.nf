process MUTRATE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(allsamples_maf)
    tuple val(meta2), path(allsamples_depths)
    tuple val(meta3), path(consensus_panel)

    output:
    tuple val(meta), path("*.mutrates.tsv"), emit: mutrates

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = task.ext.prefix ?: "${meta3.id}"
    """
    compute_mutrate.py \\
                ${allsamples_maf} \\
                ${allsamples_depths} \\
                ${consensus_panel} \\
                ${prefix}.${version};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
