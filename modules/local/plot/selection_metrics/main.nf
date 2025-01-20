process PLOT_SELECTION_METRICS {

    tag "$meta.id"
    label 'process_low'

    container 'docker.io/ferriolcalvet/deepcsa_python:20240724_latest'

    input:
    tuple val(meta), path(results_files)
    path (gene_data_df)

    output:
    tuple val(meta), path("**.pdf")  , emit: plots
    path "versions.yml"              , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_prefix = task.ext.output_prefix ?: ""
    def filters = task.ext.filters ?: ""
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
