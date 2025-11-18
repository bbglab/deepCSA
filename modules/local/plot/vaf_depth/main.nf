process PLOT_VAF_DEPTH {

    tag "$meta.id"
    label 'process_single'
    label 'time_low'
    label 'process_low_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(maf), path(mutdensity), path(omega), path(oncodrivefml)
    tuple val(meta2), path(depth)

    output:
    tuple val(meta), path("*.vaf_depth_plots.pdf"), emit: plots
    path  "versions.yml"                           , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_n = task.ext.max_n ?: 10
    def maf_arg = maf ? "--maf_file ${maf}" : ""
    def mutdensity_arg = mutdensity ? "--mutdensity_file ${mutdensity}" : ""
    def depth_arg = depth ? "--depth_file ${depth}" : ""
    def omega_arg = omega ? "--omega_file ${omega}" : ""
    def ofml_arg = oncodrivefml ? "--oncodrivefml_file ${oncodrivefml}" : ""
    """
    plot_vaf_depth.py \\
                --sample_name ${meta.id} \\
                ${maf_arg} \\
                ${mutdensity_arg} \\
                ${depth_arg} \\
                ${omega_arg} \\
                ${ofml_arg} \\
                --output_prefix ${prefix} \\
                --max_n ${max_n}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vaf_depth_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
