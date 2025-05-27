process PLOT_SELECTION_METRICS {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(results_files)
    path (gene_data_df)

    output:
    tuple val(meta), path("**.pdf")  , emit: plots
    path "versions.yml"              , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def output_prefix = task.ext.output_prefix ?: ""
    def requested_plots = task.ext.plots ?: ""
    """
    plot_selectionsideplots.py \\
                    ${prefix} \\
                    ${requested_plots} \\
                    ${gene_data_df} \\
                    ${prefix}${output_prefix} \\
                    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def output_prefix = task.ext.output_prefix ?: ""
    """
    touch ${prefix}${output_prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


//     cat > mutations_subset.conf << EOF
//     {
//         ${filters}
//     }
//     EOF

//     cat > requested_plots.conf << EOF
//     {
//         ${requested_plots}
//     }
//     EOF
