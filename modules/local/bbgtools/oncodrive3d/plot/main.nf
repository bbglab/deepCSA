process ONCODRIVE3D_PLOT {
    tag "$meta.id"
    label 'process_medium'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO pending to push the container somewhere and be able to retrieve it
    container 'docker.io/ferriolcalvet/oncodrive3d:latest'


    input:
    tuple val(meta), path(genes_csv), path(pos_csv), path(mutations_csv), path(miss_prob_json), path(seq_df_tsv)
    path(datasets)
    path(annotations)

    output:
    tuple val(meta), path("**.summary_plot.png")                               , emit: summary_plot, optional: true
    tuple val(meta), path("**.genes_plots/**.png")                             , emit: genes_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.logodds_plot.png")         , emit: logodds_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.volcano_plot.png")         , emit: volcano_plot, optional: true
    tuple val(meta), path("**.associations_plots/**.volcano_plot_gene.png")    , emit: volcano_plot_gene, optional: true
    tuple val(meta), path("**.3d_clustering_pos.annotated.csv")                , emit: pos_annotated_csv, optional: true
    tuple val(meta), path("**.log")                                            , emit: log
    path "versions.yml"                                                        , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    oncodrive3D plot -o $prefix \\
                     -g $genes_csv \\
                     -p $pos_csv \\
                     -i $mutations_csv \\
                     -m $miss_prob_json \\
                     -s $seq_df_tsv \\
                     -d $datasets \\
                     -a $annotations \\
                     -c $prefix \\
                     --output_csv \\
                     --title $prefix


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 2.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 2.0
    END_VERSIONS
    """
}

