process ONCODRIVE3D_PLOT {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/bbglab/oncodrive3d:1.0.5'

    input:
    tuple val(meta), path(genes_csv), path(pos_csv), path(mutations_csv), path(miss_prob_json), path(seq_df_tsv)
    path datasets
    path annotations

    output:
    tuple val(meta), path("**.summary_plot.png"), emit: summary_plot, optional: true
    tuple val(meta), path("**.genes_plots/**.png"), emit: genes_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.logodds_plot.png"), emit: logodds_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.volcano_plot.png"), emit: volcano_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.volcano_plot_gene.png"), emit: volcano_plot_gene, optional: true
    tuple val(meta), path("**.3d_clustering_pos.annotated.csv"), emit: pos_annotated_csv, optional: true
    tuple val(meta), path("**plot_*.log"), emit: log
    path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    oncodrive3D plot -o ${prefix} \\
                        -g ${genes_csv} \\
                        -p ${pos_csv} \\
                        -i ${mutations_csv} \\
                        -m ${miss_prob_json} \\
                        -s ${seq_df_tsv} \\
                        -d ${datasets} \\
                        -a ${annotations} \\
                        -c ${prefix} \\
                        --output_csv \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: \$(oncodrive3d --version | rev | cut -d ' ' -f 1 | rev )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: \$(oncodrive3d --version | rev | cut -d ' ' -f 1 | rev )
    END_VERSIONS
    """
}
