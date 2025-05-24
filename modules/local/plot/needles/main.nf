process PLOT_NEEDLES {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mut_files)
    path (gene_data_df)

    output:
    tuple val(meta), path("**.pdf")  , emit: plots, optional : true
    path "versions.yml"              , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_prefix = task.ext.output_prefix ?: ""
    """
    plot_needles.py \\
                    ${prefix} \\
                    ${mut_files} \\
                    ${gene_data_df} \\
                    ${prefix}${output_prefix} \\
                    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
