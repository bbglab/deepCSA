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
    """
    mkdir ${prefix}.plots
    plot_selectionsideplots.py \\
                    --sample_name ${prefix} \\
                    --outdir ${prefix}.plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

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
