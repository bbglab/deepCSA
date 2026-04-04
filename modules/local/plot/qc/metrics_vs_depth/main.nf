process PLOT_METRICS_VS_DEPTH_QC {

    tag "${group_name}"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path (all_mutdensities)
    path (depth_gene_sample)
    path (groups_json)
    val (group_name)
    path (all_adjusted_mutdensities)
    path (all_omegas_globalloc)

    output:
    path("${group_name}.metrics_depth_qc/*.pdf"), optional: true    , emit: plots
    path("${group_name}.metrics_depth_qc/*.tsv"), optional: true    , emit: tables
    path "versions.yml"                          , topic: versions

    script:
    def adjusted_arg = all_adjusted_mutdensities.name != "placeholder_no_file.tsv" ? "--adjusted-mutdensity-file ${all_adjusted_mutdensities}" : ""
    def omega_arg = all_omegas_globalloc.name != "placeholder_no_file.tsv" ? "--omegas-file ${all_omegas_globalloc}" : ""
    """
    mkdir ${group_name}.metrics_depth_qc
    metrics_vs_depth_qc.py \\
                    --mutdensity-file ${all_mutdensities} \\
                    --depth-gene-sample-file ${depth_gene_sample} \\
                    --group-definition ${groups_json} \\
                    --group-name ${group_name} \\
                    --output-dir ${group_name}.metrics_depth_qc \\
                    ${adjusted_arg} \\
                    ${omega_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${group_name}.metrics_depth_qc/${group_name}.mutdensity.depth_scatter_per_sample.pdf
    touch ${group_name}.metrics_depth_qc/${group_name}.mutdensity.depth_effect_summary.tsv
    touch ${group_name}.metrics_depth_qc/${group_name}.metrics_vs_depth_qc.status.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
