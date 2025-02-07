process PREPROCESS_DNDS {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'

    container 'docker.io/ferriolcalvet/omega:latest'

    input:
    tuple val(meta) , path(depths)
    tuple val(meta2), path (annotated_panel)


    output:
    tuple val(meta), path("*.depths_input.tsv") , emit: depths
    path "versions.yml"                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dNdS_preprocess.py ${depths} ${annotated_panel} ${prefix}.depths_input.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

