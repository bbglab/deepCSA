process PLOT_SELECTION_DEPTH {

    tag "$meta.id"
    label 'process_single'
    label 'time_low'
    label 'process_low_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(omega), path(oncodrivefml), path(depth)

    output:
    tuple val(meta), path("*.selection_depth.pdf"), emit: plots
    path  "versions.yml"                           , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def omega_arg = omega && !omega.name.equals('NO_FILE') ? "--omega_file ${omega}" : ""
    def ofml_arg = oncodrivefml && !oncodrivefml.name.equals('NO_FILE') ? "--oncodrivefml_file ${oncodrivefml}" : ""
    """
    plot_selection_depth.py \\
                --sample_name ${meta.id} \\
                ${omega_arg} \\
                ${ofml_arg} \\
                --depth_file ${depth} \\
                --output_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.selection_depth.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
