process PLOT_NEEDLES {

    tag "$meta.id"

    label 'deepcsa_core'
    label 'mem_low'

    input:
    tuple val(meta), path(mut_files)
    path (gene_data_df)

    output:
    tuple val(meta), path("**.pdf")  , emit: plots, optional : true
    path "versions.yml"              , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mkdir ${prefix}.needles ;
    plot_needles.py \\
                    --sample_name ${prefix} \\
                    --mut_file ${mut_files} \\
                    --o3d_seq_file ${gene_data_df} \\
                    --outdir ${prefix}.needles

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
