process PREPROCESS_DNDS {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta) , path(depths)
    tuple val(meta2), path (annotated_panel)


    output:
    tuple val(meta), path("*.depths_input.tsv") , emit: depths
    path "versions.yml"                         , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    dNdS_preprocess.py \\
        --depths-path ${depths} \\
        --annot-panel-path ${annotated_panel} \\
        --output ${prefix}.depths_input.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

