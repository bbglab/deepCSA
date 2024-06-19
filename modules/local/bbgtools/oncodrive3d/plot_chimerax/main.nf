process ONCODRIVE3D_PLOT_CHIMERAX {
    tag "$meta.id"
    label 'process_medium'

    // FIXME we should remove this
    queue 'bigmem'     // For some reason in the normal nodes libQt6Core.so.6 can't get imported and returns an error

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO pending to push the container somewhere and be able to retrieve it
    container 'docker.io/spellegrini87/oncodrive3d_chimerax:latest'


    input:
    tuple val(meta), path(genes_csv)
    tuple val(meta), path(pos_csv)
    tuple val(meta), path(miss_prob_json)
    tuple val(meta), path(seq_df_tsv)
    path(datasets)

    output:
    tuple val(meta), path("chimerax/attributes/**.defattr")       , emit: chimerax_defattr, optional: true
    tuple val(meta), path("chimerax/plots/**.png")                , emit: chimerax_plot, optional: true
    tuple val(meta), path("**.log")                               , emit: log
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Increase/decreaase pixel_size to decrease/increase resolution and speed/slow png generation
    """
    oncodrive3D chimerax-plot -o $prefix \\
                              -g $genes_csv \\
                              -p $pos_csv \\
                              -d $datasets \\
                              -s $seq_df_tsv \\
                              -c $prefix \\
                              --fragmented_proteins \\
                              --pixel_size 0.8 \\
                              --transparent_bg


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

