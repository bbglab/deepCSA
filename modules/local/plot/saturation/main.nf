process PLOT_SATURATION {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta) , path(results_files)  // includes all positive selection results and site comparisons
    tuple val(meta2), path(all_samples_indv_depths)
    tuple val(meta3), path(panel_file)
    path (gene_data_df)
    path (pdb_df)
    path (domains_df)
    tuple val(meta4), path (exons_depths_df)

    output:
    tuple val(meta), path("**.png")  , emit: plots
    path "versions.yml"              , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mkdir ${prefix}.plots
    plot_gene_saturation.py \\
                    --sample_name ${prefix} \\
                    --outdir ${prefix}.plots \\
                    --domain_file ${domains_df} \\
                    --exons_depths ${exons_depths_df}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
