process ONCODRIVE3D_PLOT {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO pending to push the container somewhere and be able to retrieve it
    container 'docker.io/ferriolcalvet/oncodrive3d:latest'


    input:
    tuple val(meta), path(genes_csv)
    tuple val(meta), path(pos_csv)
    tuple val(meta), path(mutations_csv)
    tuple val(meta), path(miss_prob_json)
    tuple val(meta), path(seq_df_tsv)
    path(datasets)
    path(annotations)

    output:
    tuple val(meta), path("**.summary_plot.png")                  , emit: summary_plot
    tuple val(meta), path("**.genes_plots")                       , emit: genes_plot
//    tuple val(meta), path("**.3d_clustering_pos.annotated.csv")   , emit: pos_annotated_tsv
    tuple val(meta), path("**.log")                               , emit: log
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    oncodrive3D plot -g $genes_csv \\
                     -p $pos_csv \\
                     -i $mutations_csv \\
                     -m $miss_prob_json \\
                     -s $seq_df_tsv \\
                     -d ${params.data_dir} \\
                     -a ${params.annotations_dir} \\
                     -c ${prefix} \\
                     --title ${prefix} \\
                     --c_ext \\
                     ${params.verbose ? '-v' : ''}
                    

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

