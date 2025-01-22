process PLOT_NEEDLES {

    tag "$meta.id"
    label 'process_low'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/biocontainers/seaborn:0.12.2_cv1':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seaborn:0.12.2_cv1' :
        'biocontainers/seaborn:0.12.2_cv1' }"

    input:
    tuple val(meta), path(mut_files)
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
