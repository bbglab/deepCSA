process CUSTOM_MUTATION_PROCESSING {
    tag "$meta.id"

    label 'cpu_low'
    label 'process_high_memory'
    label 'time_low'


    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta) , path(mutations_annotated)
    tuple val(meta2), path(custom_regions)

    output:
    tuple val(meta), path("*.custom.tsv")   , emit: mutations
    path "versions.yml"                     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${mutations_annotated.getBaseName()}.custom.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
