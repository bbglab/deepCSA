process CREATESAMPLEPANELS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(compact_captured_panel_annotation)
    path(depths)
    val(min_depth)

    output:
    tuple val(meta), path("*.tsv"), emit: sample_specific_panel
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    create_panel4sample.py \\
                    ${prefix}.compact.*.tsv \\
                    all_samples.depths.tsv.gz \\
                    sample_specific_panel \\
                    $min_depth

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    touch sample_specific_panel


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
