process CUSTOM_MUTATION_PROCESSING {
    tag "${meta.id}"

    label 'cpu_low'
    label 'process_high_memory'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations_annotated)
    tuple val(meta2), path(custom_regions)

    output:
    tuple val(meta), path("*.custom.tsv"), emit: mutations
    path "versions.yml", topic: versions

    script:
    // TODO reimplement python script with click
    """
    mutations_custom_processing.py \\
                    ${mutations_annotated} \\
                    ${custom_regions} \\
                    ${mutations_annotated.getBaseName()}.custom.tsv ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${mutations_annotated.getBaseName()}.custom.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
