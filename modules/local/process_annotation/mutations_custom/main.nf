process CUSTOM_MUTATION_PROCESSING {
    tag "$meta.id"

    label 'cpu_low'
    label 'process_high_memory'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta) , path(mutations_annotated)
    tuple val(meta2), path(custom_regions)

    output:
    tuple val(meta), path("*.custom.tsv")   , emit: mutations
    path "versions.yml"                     , topic: versions


    script:
    """
    mutations_custom_processing.py \\
        --mutations-file ${mutations_annotated} \\
        --custom-regions-file ${custom_regions} \\
        --output-file ${mutations_annotated.getBaseName()}.custom.tsv ;

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
