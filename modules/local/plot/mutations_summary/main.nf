process PLOT_MUTATIONS {

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

    output:
    tuple val(meta), path("*.pdf")  , emit: plots
    path "versions.yml"             , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_prefix = task.ext.output_prefix ?: ""
    def filters = task.ext.filters ?: ""
    def requested_plots = task.ext.plots ?: ""
    """
    cat > mutations_subset.conf << EOF
    {
        ${filters}
    }
    EOF

    cat > requested_plots.conf << EOF
    {
        ${requested_plots}
    }
    EOF

    plot_maf.py \\
                    ${prefix} \\
                    ${mut_files} \\
                    ${prefix}${output_prefix} \\
                    mutations_subset.conf \\
                    requested_plots.conf \\

    # plot_maf.py \\
    #                 --sample_name ${prefix} \\
    #                 --mut_file ${mut_files} \\
    #                 --out_maf ${prefix}${output_prefix}.mutations.tsv \\
    #                 --json_filters mutations_subset.conf \\
    #                 --req_fields output_formats.conf \\
    #                 ${args}
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
