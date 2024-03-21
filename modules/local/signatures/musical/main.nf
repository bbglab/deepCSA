process MUSICAL {
    tag "$meta.id"
    label 'process_high'
    label 'cpu_high'

    container 'docker.io/ferriolcalvet/musical:latest'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("**.pdf")     , emit: plots
    tuple val(meta), path("**.tsv")     , emit: stats
    tuple val(meta), path("**.pkl")     , emit: model
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_processes = task.ext.min_processes ?: "1"
    def max_processes = task.ext.max_processes ?: "10"
    """
    signatures_musical.py \\
                --sample ${prefix} \\
                --matrixfile ${matrix} \\
                --minprocess ${min_processes} \\
                --maxprocess ${max_processes} \\
                --cpus ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
