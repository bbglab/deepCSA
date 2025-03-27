process PLOT_OMEGA {

    tag "$meta.id"
    label 'process_low'

    container 'docker.io/bbglab/deepcsa-core:0.0.1-alpha'

    input:
    tuple val(meta), path(mutations), path(omegas)

    output:
    tuple val(meta), path("**.pdf")  , optional: true   , emit: plots
    path "versions.yml"                                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_prefix = task.ext.output_prefix ?: ""
    def filters = task.ext.filters ?: ""
    def requested_plots = task.ext.plots ?: "truncating,missense"
    """
    plot_selection_omega.py \\
                    ${prefix} \\
                    ${mutations} \\
                    ${omegas} \\
                    ${requested_plots} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def output_prefix = task.ext.output_prefix ?: ""
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
